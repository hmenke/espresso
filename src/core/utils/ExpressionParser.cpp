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
#include <map>
#include <string>

#include "ExpressionParser.hpp"
#include "utils/make_unique.hpp"

#include "matheval.hpp"

namespace Utils
{

class ExpressionParser::impl
{
    matheval::Parser<double> parser;
public:
    void parse(std::string const &expr)
    {
        parser.parse(expr);
    }

    void optimize()
    {
        parser.optimize();
    }

    double evaluate(std::map<std::string,double> const &st)
    {
        return parser.evaluate(st);
    }
};

ExpressionParser::ExpressionParser()
    : pimpl{Utils::make_unique<ExpressionParser::impl>()}
{}

ExpressionParser::~ExpressionParser() {}

void ExpressionParser::parse(std::string const &expr)
{
    pimpl->parse(expr);
}

void ExpressionParser::optimize()
{
    pimpl->optimize();
}

double ExpressionParser::evaluate(std::map<std::string,double> const &st)
{
    return pimpl->evaluate(st);
}

} // namespace Utils
