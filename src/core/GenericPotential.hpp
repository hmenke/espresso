#ifndef CORE_GENERIC_POTENTIAL_HPP
#define CORE_GENERIC_POTENTIAL_HPP

#include "utils/ExpressionParser.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>

#include <cassert>
#include <memory>
#include <string>

struct GenericPotential {
  double maxval = -1.0;
  std::string force_expr;
  std::string energy_expr;
  std::shared_ptr<Utils::ExpressionParser> force_parser{nullptr};
  std::shared_ptr<Utils::ExpressionParser> energy_parser{nullptr};
  bool is_parsed{false};

  void parse() {
      if (!force_parser || !energy_parser) {
          throw std::runtime_error("nullptr dereference");
      }

      force_parser->parse(force_expr);
      energy_parser->parse(energy_expr);

      is_parsed = true;
  }

  double force(double x) const {
    assert(x <= maxval);
    assert(is_parsed);
    return force_parser->evaluate({std::make_pair("r",x)});
  }

  double energy(double x) const {
    assert(x <= maxval);
    assert(is_parsed);
    return energy_parser->evaluate({std::make_pair("r",x)});
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
