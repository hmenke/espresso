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

#include "grid.hpp"
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


/** Bonded potentials from mathematical expression:
    The pair potential and pair force are given as mathematical
    expressions which are parsed and evaluated.

    @param bond_type bond type for which the interaction is defined
    @param type      type of bonded interaction: GEN_BOND_LENGTH, GEN_BOND_ANGLE
    @param max       cutoff distance
    @param energy    pair potential expression
    @param force     pair force expression

    @return 0 on success
*/
int generic_bonded_set_params(int bond_type,
                              GenericBondedInteraction type, double max,
                              std::string const &energy,
                              std::string const &force);

/* BONDED INTERACTIONS */

/** Calculate the bond length force from mathematical expression.

    @param[in]  p1       first particle
    @param[in]  p2       second particle
    @param[in]  iaparams interaction parameters
    @param[in]  dx       distance vector
    @param[out] force    bond force

    @return 0 on success, 1 on failure
*/
inline int calc_gen_bond_force(Particle *p1, Particle *p2,
                               Bonded_ia_parameters const *iaparams,
                               double dx[3], double force[3]) {
    auto const *pot = iaparams->p.gen.pot;
    auto const dist = sqrt(sqrlen(dx));

    if (dist < pot->cutoff()) {
        auto const fac = pot->force(dist) / dist;

        for (int j = 0; j < 3; j++)
            force[j] -= fac * dx[j];

        return 0;
    } else {
        return 1;
    }
}

/** Calculate the bond length energy from mathematical expression.

    @param[in]  p1       first particle
    @param[in]  p2       second particle
    @param[in]  iaparams interaction parameters
    @param[in]  dx       distance vector
    @param[out] energy   bond energy

    @return 0 on success, 1 on failure
*/
inline int gen_bond_energy(Particle *p1, Particle *p2,
                           Bonded_ia_parameters *iaparams, double dx[3],
                           double *_energy) {
    auto const *pot = iaparams->p.gen.pot;
    double dist = sqrt(sqrlen(dx));

    if (dist < pot->cutoff()) {
        *_energy = pot->energy(dist);
        return 0;
    } else {
        return 1;
    }
}


/** Calculate the bond angle bending force from mathematical expression.

    @param[in]  p_mid    middle particle
    @param[in]  p_left   left particle
    @param[in]  p_right  right particle
    @param[in]  iaparams interaction parameters
    @param[out] force1   bend force between middle and left particle
    @param[out] force2   bend force between middle and right particle

    @return 0 on success
*/
inline int calc_gen_angle_force(Particle *p_mid, Particle *p_left,
                                Particle *p_right,
                                Bonded_ia_parameters *iaparams,
                                double force1[3], double force2[3]) {
  /* vector from p_left to p_mid */
  double vec1[3];
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  double d1i = 1.0 / sqrt(sqrlen(vec1));
  for (int j = 0; j < 3; ++j)
    vec1[j] *= d1i;

  /* vector from p_mid to p_right */
  double vec2[3];
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  double d2i = 1.0 / sqrt(sqrlen(vec2));
  for (int j = 0; j < 3; ++j)
    vec2[j] *= d2i;

  /* scalar product of vec1 and vec2 */
  double cosine = scalar(vec1, vec2);
#ifdef TABANGLEMINUS // TODO: Feature guard
  double phi = acos(-cosine);
#else
  double phi = acos(cosine);
#endif
  double invsinphi = sin(phi);
  if (invsinphi < TINY_SIN_VALUE)
    invsinphi = TINY_SIN_VALUE;
  invsinphi = 1.0 / invsinphi;

  /* look up force factor */
  auto const *pot = iaparams->p.gen.pot;
  double fac = pot->force(phi);

  /* apply bend forces */
  for (int j = 0; j < 3; ++j) {
    double f1 = fac * (cosine * vec1[j] - vec2[j]) * invsinphi * d1i;
    double f2 = fac * (cosine * vec2[j] - vec1[j]) * invsinphi * d2i;
    force1[j] = (f1 - f2);
    force2[j] = -f1;
  }

  return 0;
}

/** Calculate the bond angle force on each particle from mathematical expression.

    @param[in]  p_mid    middle particle
    @param[in]  p_left   left particle
    @param[in]  p_right  right particle
    @param[in]  iaparams interaction parameters
    @param[out] force1   force on middle particle
    @param[out] force2   force on left particle
    @param[out] force3   force on right particle

    @return 0 on success
*/
inline void calc_angle_3body_generic_forces(Particle *p_mid, Particle *p_left,
                                            Particle *p_right,
                                            Bonded_ia_parameters *iaparams,
                                            double force1[3],
                                            double force2[3],
                                            double force3[3]) {
  double vec12[3]; // espresso convention
  get_mi_vector(vec12, p_mid->r.p, p_left->r.p);

  double vec21[3];
  for (int j = 0; j < 3; j++)
    vec21[j] = -vec12[j];

  double vec31[3];
  get_mi_vector(vec31, p_right->r.p, p_mid->r.p);

  double vec21_sqr = sqrlen(vec21);
  double vec21_magn = sqrt(vec21_sqr);
  double vec31_sqr = sqrlen(vec31);
  double vec31_magn = sqrt(vec31_sqr);

  double cos_phi = scalar(vec21, vec31) / (vec21_magn * vec31_magn);
  double sin_phi = sqrt(1.0 - Utils::sqr(cos_phi));

  if (cos_phi < -1.0)
    cos_phi = -TINY_COS_VALUE;
  if (cos_phi > 1.0)
    cos_phi = TINY_COS_VALUE;
#ifdef TABANGLEMINUS // TODO: Feature guard
  double phi = acos(-cos_phi);
#else
  double phi = acos(cos_phi);
#endif

  auto const *pot = iaparams->p.gen.pot;
  double dU = pot->force(phi);

  // potential dependent term (dU/dphi * 1 / sin(phi))
  double pot_dep = dU / sin_phi;

  double fj[3], fk[3];
  for (int j = 0; j < 3; ++j) {
    fj[j] =
        vec31[j] / (vec21_magn * vec31_magn) - cos_phi * vec21[j] / vec21_sqr;
    fk[j] =
        vec21[j] / (vec21_magn * vec31_magn) - cos_phi * vec31[j] / vec31_sqr;
  }

  // note that F1 = -(F2 + F3) in analytical case
  for (int j = 0; j < 3; ++j) {
    force1[j] = force1[j] - pot_dep * (fj[j] + fk[j]);
    force2[j] = force2[j] + pot_dep * fj[j];
    force3[j] = force3[j] + pot_dep * fk[j];
  }
}

/** Calculate the bond angle energy from mathematical expression.

    @param[in]  p_mid    middle particle
    @param[in]  p_left   left particle
    @param[in]  p_right  right particle
    @param[in]  iaparams interaction parameters
    @param[out] energy   bond angle energy

    @return 0 on success
*/
inline int gen_angle_energy(Particle *p_mid, Particle *p_left,
                            Particle *p_right, Bonded_ia_parameters *iaparams,
                            double *energy) {
  /* vector from p_mid to p_left */
  double vec1[3];
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  double vl1 = sqrt(sqrlen(vec1));

  /* vector from p_right to p_mid */
  double vec2[3];
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  double vl2 = sqrt(sqrlen(vec2));
/* calculate phi */
#ifdef TABANGLEMINUS // TODO: Feature guard
  double phi = acos(-scalar(vec1, vec2) / (vl1 * vl2));
#else
  double phi = acos(scalar(vec1, vec2) / (vl1 * vl2));
#endif

  auto const *pot = iaparams->p.gen.pot;
  *energy = pot->energy(phi);

  return 0;
}

#endif

#endif
