#include <gsl/gsl_errno.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/CollinearFlux.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Limits.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  const gsl_integration_fixed_type* INTEGR_TYPE = gsl_integration_fixed_jacobi;

  double unintegrated_flux(double kt2, void* params) {
    const auto args = (FluxArguments*)params;
    if (args->flux_type == KTFlux::HI_Photon_Elastic) {
      if (!args->heavy_ion)
        throw CG_FATAL("CollinearFlux") << "Heavy ion not specified!";
      return ktFlux(args->flux_type, args->x, kt2, *args->heavy_ion);
    }
    return ktFlux(args->flux_type, args->x, kt2, *args->form_factors, args->mi2, args->mf2);
  }

  CollinearFlux::CollinearFlux(const KTFlux& flux_type, const Limits& range, formfac::Parameterisation* form_fac)
      : workspace_(gsl_integration_fixed_alloc(INTEGR_TYPE, 50, range.min(), range.max(), 0., 0.)),
        params_(new FluxArguments{0., std::pow(PDG::get().mass(PDG::proton), 2), 0., flux_type, form_fac, nullptr}),
        function_({&unintegrated_flux, (void*)params_.get()}) {}

  CollinearFlux::CollinearFlux(const KTFlux& flux_type, const Limits& range, HeavyIon* hi)
      : workspace_(gsl_integration_fixed_alloc(INTEGR_TYPE, 50, range.min(), range.max(), 0., 0.)),
        params_(new FluxArguments{0., std::pow(PDG::get().mass(PDG::proton), 2), 0., flux_type, nullptr, hi}),
        function_({&unintegrated_flux, (void*)params_.get()}) {}

  double CollinearFlux::operator()(double x, double mx) const {
    double result = 0.;
    params_->x = x;
    params_->mf2 = mx * mx;
    const int res = gsl_integration_fixed(&function_, &result, workspace_.get());
    if (res != GSL_SUCCESS)
      CG_ERROR("CollinearFlux") << gsl_strerror(res);
    result *= M_1_PI;
    return result;
  }
}  // namespace cepgen
