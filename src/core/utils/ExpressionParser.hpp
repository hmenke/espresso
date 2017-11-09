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
#ifndef CORE_UTILS_EXPRESSION_PARSER_HPP
#define CORE_UTILS_EXPRESSION_PARSER_HPP

#include <map>
#include <memory>
#include <string>

namespace Utils
{

/** @brief Parse a mathematical expression
 *
 * This can parse and evaluate a mathematical expression for a given
 * symbol table using Boost Matheval.  The templates of Boost.Spirit
 * (which is the basis of Boost Matheval) are very expensive to parse
 * and instantiate, which is why we hide it behind an opaque pointer.
 *
 * The drawback of this approach is that calls can no longer be
 * inlined and because the pointer crosses translation unit
 * boundaries, dereferencing it can also not be optimized out at
 * compile time.  We have to rely entirely on link-time optimization
 * which might be not as good.
 *
 * The pointer to the implementation is a std::unique_ptr which makes
 * the class not copyable but only moveable.  Copying shouldn't be
 * required but is easy to implement.
 */
class ExpressionParser
{
    class impl;
    std::unique_ptr<impl> pimpl;
public:
    /** @brief Constructor */
    ExpressionParser();

    /** @brief Destructor */
    ~ExpressionParser();

    /** @brief Parse the mathematical expression into an abstract syntax tree
     *
     * @param[in] expr The expression given as a std::string
     */
    void parse(std::string const &expr);

    /** @brief Evaluate the abstract syntax tree for a given symbol table
     *
     * @param[in] st The symbol table
     */
    double evaluate(std::map<std::string,double> const &st);
};

} // namespace Utils

#endif // CORE_UTILS_EXPRESSION_PARSER_HPP
