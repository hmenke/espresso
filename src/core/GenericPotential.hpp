#ifndef CORE_GENERIC_POTENTIAL_HPP
#define CORE_GENERIC_POTENTIAL_HPP

#include "utils/ExpressionParser.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>

#include <cassert>
#include <string>

struct GenericPotential {
  double maxval = -1.0;
  std::string force_expr;
  std::string energy_expr;

  double force(double x) const {
    assert(x <= maxval);
    Utils::ExpressionParser parser;
    parser.parse(force_expr);
    return parser.evaluate({std::make_pair("r",x)});
  }

  double energy(double x) const {
    assert(x <= maxval);
    Utils::ExpressionParser parser;
    parser.parse(energy_expr);
    return parser.evaluate({std::make_pair("r",x)});
  }

  double cutoff() const { return maxval; }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &maxval;
    ar &force_expr;
    ar &energy_expr;
  }
};

#endif
