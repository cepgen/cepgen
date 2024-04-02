/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2024  Laurent Forthomme
 *                     2016  Antoni Szczurek
 *                     2016  Volodymyr Uleshchenko
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

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Message.h"

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
      explicit SzczurekUleshchenko(const ParametersList& params)
          : Parameterisation(params), q2_shift_(steerAs<double, float>("q2shift")) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Szczurek-Uleshchenko (based on GRV parton content)");
        desc.add<double>("q2shift", 0.8);
        return desc;
      }

    private:
      void eval() override {
        auto amu2 = (float)args_.q2 + q2_shift_;  // shift the overall scale
        float xuv, xdv, xus, xds, xss, xg;
        auto xbj_arg = (float)args_.xbj;

        grv95lo_(xbj_arg, amu2, xuv, xdv, xus, xds, xss, xg);

        CG_DEBUG_LOOP("SzczurekUleshchenko")
            << "Form factor content at xB = " << args_.xbj << " (scale = " << amu2 << " GeV^2):\n\t"
            << "  valence quarks: u / d     = " << xuv << " / " << xdv << "\n\t"
            << "  sea quarks:     u / d / s = " << xus << " / " << xds << " / " << xss << "\n\t"
            << "  gluons:                   = " << xg;

        // standard partonic structure function
        const double F2_aux = 4. / 9. * (xuv + 2. * xus) + 1. / 9. * (xdv + 2. * xds) + 1. / 9. * (2. * xss);
        setF2(F2_aux * args_.q2 / amu2);  // F2 corrected for low Q^2 behaviour
      }
      const float q2_shift_;  ///< \f$Q^2\f$ scale shift
    };
  }  // namespace strfun
}  // namespace cepgen
using cepgen::strfun::SzczurekUleshchenko;
REGISTER_STRFUN("SzczurekUleshchenko", 12, SzczurekUleshchenko);
