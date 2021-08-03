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

#include <gsl/gsl_errno.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/CollinearFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Limits.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  double unintegrated_flux(double kt2, void* params) {
    const auto args = static_cast<FluxArguments*>(params);
    if (args->flux_type == KTFlux::HI_Photon_Elastic) {
      if (!args->heavy_ion)
        throw CG_FATAL("CollinearFlux") << "Heavy ion not specified!";
      return ktFlux(args->flux_type, args->x, kt2, *args->heavy_ion) / kt2;
    }
    return ktFlux(args->flux_type, args->x, kt2, *args->form_factors, args->mi2, args->mf2) / kt2;
  }

  CollinearFlux::CollinearFlux(formfac::Parameterisation* form_fac, const Limits& kt2_range)
      : workspace_(
            gsl_integration_fixed_alloc(gsl_integration_fixed_jacobi, 50, kt2_range.min(), kt2_range.max(), 0., 0.)),
        params_(
            new FluxArguments{0., std::pow(PDG::get().mass(PDG::proton), 2), 0., KTFlux::invalid, form_fac, nullptr}),
        function_({&unintegrated_flux, (void*)params_.get()}) {}

  CollinearFlux::CollinearFlux(HeavyIon* hi, const Limits& range)
      : workspace_(gsl_integration_fixed_alloc(gsl_integration_fixed_jacobi, 50, range.min(), range.max(), 0., 0.)),
        params_(new FluxArguments{0., std::pow(PDG::get().mass(PDG::proton), 2), 0., KTFlux::invalid, nullptr, hi}),
        function_({&unintegrated_flux, (void*)params_.get()}) {}

  double CollinearFlux::operator()(double x, double mx, const KTFlux& flux) const {
    double result = 0.;
    params_->x = x;
    params_->mf2 = mx * mx;
    params_->flux_type = flux;
    const int res = gsl_integration_fixed(&function_, &result, workspace_.get());
    if (res != GSL_SUCCESS)
      CG_ERROR("CollinearFlux") << gsl_strerror(res);
    result *= M_1_PI;
    return result;
  }
}  // namespace cepgen
