/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGen_Integration_FunctionIntegrand_h
#define CepGen_Integration_FunctionIntegrand_h

#include <functional>

#include "CepGen/Integration/Integrand.h"

namespace cepgen {
  /// Wrapper to the function to be integrated
  class FunctionIntegrand : public Integrand {
  public:
    explicit FunctionIntegrand(size_t, const std::function<double(const std::vector<double>&)>&);

    double eval(const std::vector<double>&) override;
    size_t size() const override { return num_dimensions_; }

  private:
    std::function<double(const std::vector<double>&)> function_;
    size_t num_dimensions_;
  };
}  // namespace cepgen

#endif
