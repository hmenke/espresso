/*
  Copyright (C) 2017 The ESPResSo project

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

#define BOOST_TEST_MODULE Expression parser test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

#include "utils/ExpressionParser.hpp"

// We only test the integration into ESPResSo as Boost Matheval has
// plenty of tests itself.

BOOST_AUTO_TEST_CASE(integration1) {
    std::string expr = "cbrt(x/2 + sqrt(x**2/4 + y**3/24))";
    double x = 2.0, y = -1.0;

    Utils::ExpressionParser parser;
    BOOST_CHECK_NO_THROW(parser.parse(expr));

    std::map<std::string,double> symbol_table = {
        std::make_pair("x",x),
        std::make_pair("y",y),
    };

    double result = 0;
    BOOST_CHECK_NO_THROW(result = parser.evaluate(symbol_table));

    double expected = std::cbrt(x/2. + std::sqrt(std::pow(x,2.)/4. + std::pow(y,3.)/24.));
    BOOST_CHECK_CLOSE_FRACTION(result, expected, std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(integration2) {
    std::string expr = "(";

    Utils::ExpressionParser parser;

    // Parsing should fail
    BOOST_CHECK_THROW(parser.parse(expr), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(integration3) {
    std::string expr = "x";

    Utils::ExpressionParser parser;

    BOOST_CHECK_NO_THROW(parser.parse(expr));

    // Evaluating should fail
    BOOST_CHECK_THROW(parser.evaluate({}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(integration4) {
    Utils::ExpressionParser parser;

    // Missing expression evaluates to zero
    BOOST_CHECK_EQUAL(parser.evaluate({}), 0);
}
