/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file
 *  This file contains everything needed to create a start-up configuration
 *  of (partially charged) polymer chains with counterions and salt molecules,
 *  assigning velocities to the particles and cross-linking the polymers if
 *  necessary.
 *
 *  The corresponding header file is polymer.hpp.
 */

#include "polymer.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "constraints.hpp"
#include "constraints/ShapeBasedConstraint.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/vec_rotate.hpp>

#include <cmath>
#include <random>
#include <stdexcept>

template <class RNG> static Utils::Vector3d random_position(RNG &rng) {
  Utils::Vector3d v;
  for (int i = 0; i < 3; ++i)
    v[i] = box_geo.length()[i] * rng();
  return v;
}

template <class RNG> static Utils::Vector3d random_unit_vector(RNG &rng) {
  Utils::Vector3d v;
  double const phi = acos(1. - 2. * rng());
  double const theta = 2. * Utils::pi() * rng();
  v[0] = sin(phi) * cos(theta);
  v[1] = sin(phi) * sin(theta);
  v[2] = cos(phi);
  v /= v.norm();
  return v;
}

/** Determines whether a given position @p pos is valid, i.e., it doesn't
 *  collide with existing or buffered particles, nor with existing constraints
 *  (if @c respect_constraints).
 *  @param pos                   the trial position in question
 *  @param positions             buffered positions to respect
 *  @param partCfg               existing particles to respect
 *  @param min_distance          threshold for the minimum distance between
 *                               trial position and buffered/existing particles
 *  @param respect_constraints   whether to respect constraints
 *  @return true if valid position, false if not.
 */
static bool
is_valid_position(Utils::Vector3d const *pos,
                  std::vector<std::vector<Utils::Vector3d>> const *positions,
                  PartCfg &partCfg, double const min_distance,
                  int const respect_constraints) {
  Utils::Vector3d const folded_pos = folded_position(*pos, box_geo);
  // check if constraint is violated
  if (respect_constraints) {
    for (auto &c : Constraints::constraints) {
      auto cs =
          std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
      if (cs) {
        double d;
        Utils::Vector3d v;

        cs->calc_dist(folded_pos, d, v);

        if (d <= 0) {
          return false;
        }
      }
    }
  }

  if (min_distance > 0) {
    // check for collision with existing particles
    if (distto(partCfg, *pos, -1) < min_distance) {
      return false;
    }
    // check for collision with buffered positions
    double buff_mindist = std::numeric_limits<double>::infinity();
    double h;
    for (auto const &p : *positions) {
      for (auto const &m : p) {
        h = (folded_position(*pos, box_geo) - folded_position(m, box_geo))
                .norm2();
        buff_mindist = std::min(h, buff_mindist);
      }
    }
    if (std::sqrt(buff_mindist) < min_distance) {
      return false;
    }
  }
  return true;
}

std::vector<std::vector<Utils::Vector3d>>
draw_polymer_positions(PartCfg &partCfg, int const n_polymers,
                       int const beads_per_chain, double const bond_length,
                       std::vector<Utils::Vector3d> const &start_positions,
                       double const min_distance, int const max_tries,
                       int const use_bond_angle, double const bond_angle,
                       int const respect_constraints, int const seed) {
  std::vector<std::vector<Utils::Vector3d>> positions(
      n_polymers, std::vector<Utils::Vector3d>(beads_per_chain));

  auto rng = [mt = std::mt19937{static_cast<unsigned>(seed)},
              dist = std::uniform_real_distribution<double>(
                  0.0, 1.0)]() mutable { return dist(mt); };

  Utils::Vector3d trial_pos;
  int attempts_mono, attempts_poly;

  // make sure that if given, all starting positions are valid
  if ((not start_positions.empty()) and
      std::any_of(start_positions.begin(), start_positions.end(),
                  [&positions, &partCfg, min_distance,
                   respect_constraints](Utils::Vector3d v) {
                    return not is_valid_position(&v, &positions, partCfg,
                                                 min_distance,
                                                 respect_constraints);
                  }))
    throw std::runtime_error("Invalid start positions.");
  // use (if none given, random) starting positions for every first monomer
  for (int p = 0; p < n_polymers; ++p) {
    attempts_mono = 0;
    // first monomer for all polymers
    if (start_positions.empty()) {
      do {
        trial_pos = random_position(rng);
        attempts_mono++;
      } while ((not is_valid_position(&trial_pos, &positions, partCfg,
                                      min_distance, respect_constraints)) and
               (attempts_mono < max_tries));
      if (attempts_mono == max_tries) {
        throw std::runtime_error("Failed to create polymer start positions.");
      }
      positions[p][0] = trial_pos;
    } else {
      positions[p][0] = start_positions[p];
    }
  }

  // create remaining monomers' positions
  for (int p = 0; p < n_polymers; ++p) {
    attempts_poly = 0;
    for (int m = 1; m < beads_per_chain; ++m) {
      attempts_mono = 0;
      do {
        if (m == 0) {
          // m == 0 is only true after a failed attempt to position a polymer
          trial_pos = random_position(rng);
        } else if (not use_bond_angle or m < 2) {
          // random step, also necessary if angle is set, placing the second
          // bead
          trial_pos =
              positions[p][m - 1] + bond_length * random_unit_vector(rng);
        } else {
          // use prescribed angle
          Utils::Vector3d last_vec = positions[p][m - 1] - positions[p][m - 2];
          trial_pos = positions[p][m - 1] +
                      Utils::vec_rotate(
                          vector_product(last_vec, random_unit_vector(rng)),
                          bond_angle, -last_vec);
        }
        attempts_mono++;
      } while ((not is_valid_position(&trial_pos, &positions, partCfg,
                                      min_distance, respect_constraints)) and
               (attempts_mono < max_tries));

      if (attempts_mono == max_tries) {
        if (attempts_poly < max_tries) {
          // if start positions have to be respected: fail ...
          if (start_positions.empty()) {
            throw std::runtime_error("Failed to create polymer positions with "
                                     "given start positions.");
          }
          // ... otherwise retry to position the whole polymer
          attempts_poly++;
          m = -1;
        } else {
          // ... but only if max_tries has not been exceeded.
          throw std::runtime_error("Failed to create polymer positions.");
        }
      } else {
        positions[p][m] = trial_pos;
      }
    }
  }
  return positions;
}
