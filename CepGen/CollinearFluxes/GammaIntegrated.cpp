/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2022  Laurent Forthomme
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

#include "CepGen/CollinearFluxes/Parameterisation.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Integration/AnalyticIntegrator.h"
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/FunctionsWrappers.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  namespace collflux {
    struct FluxArguments {
      double x{0.}, mi2{0.}, mf2{0.};
      Beam::KTFlux flux_type{Beam::KTFlux::P_Photon_Elastic};
      formfac::Parameterisation* form_factors{nullptr};
      strfun::Parameterisation* structure_functions{nullptr};
      HeavyIon* heavy_ion{nullptr};
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
            form_fac_(formfac::FormFactorsFactory::get().build(steer<std::string>("formFactors"))),
            integr_(AnalyticIntegratorFactory::get().build(params.get<ParametersList>("analyticalIntegrator"))),
            params_(new FluxArguments{0., mp2_, 0., flux_, form_fac_.get(), nullptr, nullptr}) {
        const auto& plist_strfun = steer<ParametersList>("structureFunctions");
        if (!plist_strfun.empty()) {
          str_fun_ = strfun::StructureFunctionsFactory::get().build(plist_strfun);
          params_->structure_functions = str_fun_.get();
        }
        func_.reset(new utils::Function1D(unintegrated_flux));
        CG_INFO("GammaIntegrated").log([&](auto& log) {
          log << "kt flux-integrated collinear flux evaluator initialised.\n\t"
              << "Q^2 integration range: " << q2_range_ << " GeV^2\n\t"
              << "Nucleon/HI: " << hi_ << "\n\t";
          if (str_fun_)
            log << "Structure functions modelling: " << *str_fun_ << "\n\t";
          log << "Form factors modelling: " << *form_fac_ << ".";
        });
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("kt-integrated photon flux");
        desc.addAs<int, Beam::KTFlux>("ktFlux", Beam::KTFlux::P_Photon_Elastic);
        desc.addAs<pdgid_t, HeavyIon>("heavyIon", HeavyIon::proton());
        desc.add<std::string>("formFactors", formfac::gFFStandardDipoleHandler);
        desc.add<ParametersDescription>("structureFunctions",
                                        strfun::StructureFunctionsFactory::get().describeParameters(11));
        desc.add<ParametersDescription>("analyticalIntegrator", ParametersDescription().setName<std::string>("gsl"))
            .setDescription("Steering parameters for the analytical integrator");
        return desc;
      }

      double operator()(double x, double mx) const override {
        static const Limits x_valid_range{0., 1.};
        if (x_valid_range.contains(x))
          return 0.;
        params_->x = x;
        params_->mf2 = mx * mx;
        return 2. * M_PI * integr_->eval(*func_, params_.get(), q2_range_) / x;
      }

    private:
      const Beam::KTFlux flux_;
      const HeavyIon hi_;
      std::unique_ptr<formfac::Parameterisation> form_fac_;
      std::unique_ptr<strfun::Parameterisation> str_fun_;
      std::unique_ptr<utils::Function1D> func_;
      std::unique_ptr<AnalyticIntegrator> integr_;
      std::unique_ptr<FluxArguments> params_;
    };
  }  // namespace collflux
}  // namespace cepgen
typedef cepgen::collflux::GammaIntegrated CF_GI;
REGISTER_COLLFLUX("GammaIntegrated", CF_GI)
