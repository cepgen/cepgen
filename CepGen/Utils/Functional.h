/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Utils_Functional_h
#define CepGen_Utils_Functional_h

#include <array>
#include <string>
#include <vector>

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class ParametersList;
  namespace utils {
    /// A string-to-functional parser
    /// \author L. Forthomme <laurent.forthomme@cern.ch>
    /// \date 21 Aug 2017
    class Functional : public NamedModule<Functional, std::string> {
    public:
      /// Default constructor
      explicit Functional(const ParametersList&);

      /// Build a collection of parameters to define a functional from its mathematical expression
      /// \param[in] expr Mathematical expression to evaluate
      /// \param[in] vars List of expression variables
      static ParametersList fromExpression(const std::string& expr, const std::vector<std::string>& vars);
      static ParametersDescription description();

      /// Compute the functional for a given value of the variable (one-dimensional case)
      /// \param[in] x Variable value
      double operator()(double x) const;
      /// Compute the functional for a given value of the variables
      /// \param[in] x Variables values
      double operator()(const std::vector<double>& x) const;

      /// List of string variable names
      const std::vector<std::string>& variables() const { return vars_orig_; }
      /// String expression held by this functional parser
      const std::string& expression() const { return expression_orig_; }

    protected:
      /// Compute the functional for a given value of the variables
      virtual double eval() const = 0;

    private:
      std::vector<std::string> vars_orig_;  ///< User-defined variables to be reached
      std::string expression_orig_;         ///< User-defined expression

    protected:
      std::vector<std::string> vars_;       ///< Computer-readable variable to be reached
      std::string expression_;              ///< Computer-readable expression
      mutable std::vector<double> values_;  ///< Last arguments list fed to the functional
    };
  }  // namespace utils
}  // namespace cepgen

#endif
