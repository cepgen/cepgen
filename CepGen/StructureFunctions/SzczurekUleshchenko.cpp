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

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace {
  extern "C" {
  extern void grv95lo_(float&, float&, float&, float&, float&, float&, float&, float&);
  }
}  // namespace

namespace cepgen {
  namespace strfun {
    /// Szczurek and Uleshchenko modelling of \f$F_2\f$ based on GRV parton content \cite Szczurek:1999wp
    class SzczurekUleshchenko final : public Parameterisation {
    public:
      explicit SzczurekUleshchenko(const ParametersList&);
      SzczurekUleshchenko& eval(double xbj, double q2) override;

      static ParametersDescription description();

    private:
      /// \f$Q^2\f$ scale shift
      const float q2_shift_;
    };

    SzczurekUleshchenko::SzczurekUleshchenko(const ParametersList& params)
        : Parameterisation(params), q2_shift_(params.getAs<double, float>("q2shift")) {}

    SzczurekUleshchenko& SzczurekUleshchenko::eval(double xbj, double q2) {
      auto amu2 = (float)q2 + q2_shift_;  // shift the overall scale
      float xuv, xdv, xus, xds, xss, xg;
      auto xbj_arg = (float)xbj;

      grv95lo_(xbj_arg, amu2, xuv, xdv, xus, xds, xss, xg);

      CG_DEBUG_LOOP("SzczurekUleshchenko")
          << "Form factor content at xB = " << xbj << " (scale = " << amu2 << " GeV^2):\n\t"
          << "  valence quarks: u / d     = " << xuv << " / " << xdv << "\n\t"
          << "  sea quarks:     u / d / s = " << xus << " / " << xds << " / " << xss << "\n\t"
          << "  gluons:                   = " << xg;

      // standard partonic structure function
      const double F2_aux = 4. / 9. * (xuv + 2. * xus) + 1. / 9. * (xdv + 2. * xds) + 1. / 9. * (2. * xss);

      F2 = F2_aux * q2 / amu2;  // F2 corrected for low Q^2 behaviour

      return *this;
    }

    ParametersDescription SzczurekUleshchenko::description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Szczurek-Uleshchenko modelling of F2 based on GRV parton content");
      desc.add<double>("q2shift", 0.8);
      return desc;
    }
  }  // namespace strfun
}  // namespace cepgen

REGISTER_STRFUN(strfun::Type::SzczurekUleshchenko, SzczurekUleshchenko, strfun::SzczurekUleshchenko)
