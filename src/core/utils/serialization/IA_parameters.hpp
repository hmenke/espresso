#ifndef SERIALIZATION_IA_PARAMETERS_HPP
#define SERIALIZATION_IA_PARAMETERS_HPP

#include <memory>

#include "utils/ExpressionParser.hpp"
#include "core/interaction_data.hpp"

namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar, IA_parameters &p,
          const unsigned int /* file_version */) {
  ar >> make_array(reinterpret_cast<char *>(&p), sizeof(IA_parameters));

#ifdef TABULATED
  TabulatedPotential tab;
  ar >> tab;

  new (&(p.TAB)) TabulatedPotential(std::move(tab));
#endif

#if true // TODO: Feature guard
  GenericPotential gen;
  ar >> gen;

  gen.force_parser = std::make_shared<Utils::ExpressionParser>();
  gen.energy_parser = std::make_shared<Utils::ExpressionParser>();

  gen.parse();

  new (&(p.GEN)) GenericPotential(std::move(gen));
#endif
}

template <typename Archive>
void save(Archive &ar, IA_parameters const &p,
          const unsigned int /* file_version */) {
  ar << make_array(reinterpret_cast<char const *>(&p), sizeof(IA_parameters));

#ifdef TABULATED
  ar << p.TAB;
#endif

#if true // TODO: Feature guard
  ar << p.GEN;
#endif
}

template <class Archive>
void serialize(Archive &ar, IA_parameters &p, const unsigned int file_version) {
  split_free(ar, p, file_version);
}
}
}

#endif
