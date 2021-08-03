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

#include <cmath>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

namespace cepgen {
  namespace formfac {
    class StandardDipole final : public Parameterisation {
    public:
      using Parameterisation::Parameterisation;
      static std::string description() { return "Standard dipole"; }

    private:
      void compute(double q2) override {
        GE = pow(1. + q2 * inv_sq_scale_param_, -2.);
        GM = MU * GE;
      }

      static constexpr double inv_sq_scale_param_ = 1. / 0.71;  // for r_p = 0.81 fm
      //static constexpr double inv_sq_scale_param_ = 1. / 0.66; // for r_p = 0.84 fm
    };
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FF_MODEL(gFFStandardDipoleHandler, StandardDipole)
