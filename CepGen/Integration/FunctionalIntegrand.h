/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#ifndef CepGen_Integration_FunctionalIntegrand_h
#define CepGen_Integration_FunctionalIntegrand_h

#include <memory>

#include "CepGen/Integration/Integrand.h"
#include "CepGen/Utils/Functional.h"

namespace cepgen {
  /// Wrapper to a string-built functional to be integrated
  class FunctionalIntegrand final : public Integrand {
  public:
    explicit FunctionalIntegrand(const std::string& expression,
                                 const std::vector<std::string>& variables,
                                 const std::string& functional_evaluator);

    double eval(const std::vector<double>&) override;
    size_t size() const override;

  private:
    std::unique_ptr<utils::Functional> functional_;
  };
}  // namespace cepgen

#endif
