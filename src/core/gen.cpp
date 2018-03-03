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
#include "gen.hpp"

#ifdef EXPRESSION

#include <memory>

#include "communication.hpp"
#include "utils/ExpressionParser.hpp"

int generic_set_params(int part_type_a, int part_type_b,
                       double max, std::string const &energy,
                       std::string const &force) {
  auto data = get_ia_param_safe(part_type_a, part_type_b);

  data->GEN.maxval = max;

  data->GEN.force_expr = force;
  data->GEN.energy_expr = energy;

  data->GEN.force_parser = std::make_shared<Utils::ExpressionParser>();
  data->GEN.energy_parser = std::make_shared<Utils::ExpressionParser>();

  data->GEN.parse();

  mpi_bcast_ia_params(part_type_a, part_type_b);

  return 0;
}


int generic_bonded_set_params(int bond_type,
                              GenericBondedInteraction type, double max,
                              std::string const &energy,
                              std::string const &force) {
  if (bond_type < 0)
    throw std::runtime_error("Ill-formed bond type.");

  make_bond_type_exist(bond_type);

  /* set types */
  bonded_ia_params[bond_type].type = BONDED_IA_GENERIC;
  bonded_ia_params[bond_type].p.gen.type = type;
  bonded_ia_params[bond_type].p.gen.pot = new GenericPotential;
  auto pot = bonded_ia_params[bond_type].p.gen.pot;

  /* set number of interaction partners */
  if (type == GEN_BOND_LENGTH) {
    pot->maxval = max;
    bonded_ia_params[bond_type].num = 1;
  } else if (type == GEN_BOND_ANGLE) {
    pot->maxval = PI + ROUND_ERROR_PREC;
    bonded_ia_params[bond_type].num = 2;
  } else {
    throw std::runtime_error("Unsupported generic bond type.");
  }

  pot->force_expr = force;
  pot->energy_expr = energy;

  pot->force_parser = std::make_shared<Utils::ExpressionParser>();
  pot->energy_parser = std::make_shared<Utils::ExpressionParser>();

  pot->parse();

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

#endif
