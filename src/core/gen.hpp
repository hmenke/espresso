/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CORE_GENERIC_HPP
#define CORE_GENERIC_HPP

#include "config.hpp"

#if true // TODO: Feature guard

#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

/** Non-Bonded potentials from mathematical expression:
    The pair potential and pair force are given as mathematical
    expressions which are parsed and evaluated.

    @param part_type_a particle type for which the interaction is defined
    @param part_type_b particle type for which the interaction is defined
    @param max         cutoff distance
    @param energy      pair potential expression
    @param force       pair force expression

    @return 0 on success
*/
int generic_set_params(int part_type_a, int part_type_b,
                       double max, std::string const &energy,
                       std::string const &force);

/** Add a non-bonded pair force from a mathematical expression. */
inline void add_generic_pair_force(const Particle *const p1,
                                   const Particle *const p2,
                                   IA_parameters const *ia_params,
                                   double const d[3], double dist,
                                   double force[3]) {
    if (dist < ia_params->GEN.cutoff()) {
        auto const fac = ia_params->GEN.force(dist) / dist;

        for (int j = 0; j < 3; j++)
            force[j] -= fac * d[j];
    }
}

/** Add a non-bonded pair energy from a mathematical expression. */
inline double generic_pair_energy(Particle const *, Particle const *,
                                  IA_parameters const *ia_params,
                                  const double d[3], double dist) {
    if (dist < ia_params->GEN.cutoff()) {
        return ia_params->GEN.energy(dist);
    } else {
        return 0.0;
    }
}

#endif

#endif
