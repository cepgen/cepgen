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

#ifndef CepGen_Physics_CollinearFlux_h
#define CepGen_Physics_CollinearFlux_h

#include <gsl/gsl_integration.h>

#include <memory>

#include "CepGen/Physics/Beam.h"

namespace cepgen {
  class Limits;
  namespace formfac {
    class Parameterisation;
  }
  namespace strfun {
    class Parameterisation;
  }
  class HeavyIon;
  struct FluxArguments {
    double x, mi2, mf2;
    Beam::KTFlux flux_type;
    formfac::Parameterisation* form_factors;
    strfun::Parameterisation* structure_functions;
    HeavyIon* heavy_ion;
  };
  class CollinearFlux {
  public:
    explicit CollinearFlux(formfac::Parameterisation*, strfun::Parameterisation*, const Limits&);
    explicit CollinearFlux(HeavyIon* hi, const Limits& kt2_range);

    double operator()(double x, double mx, const Beam::KTFlux& flux_type) const;

  private:
    struct gsl_integration_fixed_workspace_del {
      void operator()(gsl_integration_fixed_workspace* int_wsp) { gsl_integration_fixed_free(int_wsp); }
    };
    std::unique_ptr<gsl_integration_fixed_workspace, gsl_integration_fixed_workspace_del> workspace_;
    std::unique_ptr<FluxArguments> params_;
    mutable gsl_function function_;
  };
}  // namespace cepgen

#endif
