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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace collflux {
    struct FluxArguments {
      double x, mi2, mf2;
      Beam::KTFlux flux_type;
      formfac::Parameterisation* form_factors;
      strfun::Parameterisation* structure_functions;
      HeavyIon* heavy_ion;
    };

    double unintegrated_flux(double kt2, void* params) {
      const auto args = static_cast<FluxArguments*>(params);
      if (args->flux_type == Beam::KTFlux::HI_Photon_Elastic) {
        if (!args->heavy_ion)
          throw CG_FATAL("CollinearFlux") << "Heavy ion not specified!";
        return Beam::ktFluxHI(args->flux_type, args->x, kt2, *args->heavy_ion) / kt2;
      }
      return Beam::ktFluxNucl(
                 args->flux_type, args->x, kt2, args->form_factors, args->structure_functions, args->mi2, args->mf2) /
             kt2;
    }

    class GammaIntegrated : public Parameterisation {
    public:
      explicit GammaIntegrated(const ParametersList& params)
          : Parameterisation(params),
            flux_(steerAs<int, Beam::KTFlux>("ktFlux")),
            hi_(steerAs<pdgid_t, HeavyIon>("heavyIon")),
            form_fac_(formfac::FormFactorsFactory::get().build(steer<ParametersList>("formFactors"))),
            str_fun_(strfun::StructureFunctionsFactory::get().build(steer<ParametersList>("structureFunctions"))),
            workspace_(
                gsl_integration_fixed_alloc(gsl_integration_fixed_jacobi, 50, t_range_.min(), t_range_.max(), 0., 0.)),
            params_(new FluxArguments{0., mp2_, 0., flux_, form_fac_.get(), str_fun_.get(), nullptr}),
            function_({&unintegrated_flux, (void*)params_.get()}) {
        CG_INFO("GammaIntegrated") << "kt flux-integrated collinear flux evaluator initialised.\n\t"
                                   << "Q^2 integration range: " << t_range_ << " GeV^2\n\t"
                                   << "Nucleon/HI: " << hi_ << "\n\t"
                                   << "Structure functions modelling: " << *str_fun_ << "\n\t"
                                   << "Form factors modelling: " << *form_fac_ << ".";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("kt-integrated photon flux");
        desc.addAs<int, Beam::KTFlux>("ktFlux", Beam::KTFlux::P_Photon_Elastic);
        desc.add<pdgid_t>("heavyIon", HeavyIon::proton());
        desc.add<ParametersDescription>("formFactors", ParametersDescription());
        desc.add<ParametersDescription>("structureFunctions", ParametersDescription());
        return desc;
      }

      double operator()(double x, double mx) const override {
        double result = 0.;
        params_->x = x;
        params_->mf2 = mx * mx;
        const int res = gsl_integration_fixed(&function_, &result, workspace_.get());
        if (res != GSL_SUCCESS)
          CG_ERROR("CollinearFlux") << gsl_strerror(res);
        result *= M_1_PI;
        return result;
      }

    private:
      const Beam::KTFlux flux_;
      const HeavyIon hi_{HeavyIon::proton()};
      std::unique_ptr<formfac::Parameterisation> form_fac_;
      std::unique_ptr<strfun::Parameterisation> str_fun_;
      struct gsl_integration_fixed_workspace_del {
        void operator()(gsl_integration_fixed_workspace* int_wsp) { gsl_integration_fixed_free(int_wsp); }
      };
      std::unique_ptr<gsl_integration_fixed_workspace, gsl_integration_fixed_workspace_del> workspace_;
      std::unique_ptr<FluxArguments> params_;
      mutable gsl_function function_;
    };
  }  // namespace collflux
}  // namespace cepgen

REGISTER_COLLFLUX(GammaIntegrated, collflux::GammaIntegrated);
