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
    class GammaIntegrated : public Parameterisation {
    public:
      explicit GammaIntegrated(const ParametersList& params)
          : Parameterisation(params),
            flux_(steerAs<int, Beam::KTFlux>("ktFlux")),
            hi_(steerAs<pdgid_t, HeavyIon>("heavyIon")),
            form_fac_(formfac::FormFactorsFactory::get().build(steer<std::string>("formFactors"))),
            integr_(AnalyticIntegratorFactory::get().build(params.get<ParametersList>("analyticalIntegrator"))) {
        const auto& plist_strfun = steer<ParametersList>("structureFunctions");
        if (!plist_strfun.empty())
          str_fun_ = strfun::StructureFunctionsFactory::get().build(plist_strfun);

        // initialise the function to integrate
        if (hi_ != HeavyIon::proton() && hi_ != HeavyIon::neutron())
          func_.reset(new utils::Function1D([&](double kt2, void* params) {
            if (kt2 < 0.)
              return 0.;
            if (!params)
              throw CG_FATAL("CollinearFlux") << "Invalid parameters block fed to the integrand!";
            const auto& args = *static_cast<FluxArguments*>(params);
            return Beam::ktFluxHI(flux_, args.x, kt2, hi_) / kt2;
          }));
        else
          func_.reset(new utils::Function1D([&](double kt2, void* params) {
            if (kt2 < 0.)
              return 0.;
            if (!params)
              throw CG_FATAL("CollinearFlux") << "Invalid parameters block fed to the integrand!";
            const auto& args = *static_cast<FluxArguments*>(params);
            return Beam::ktFluxNucl(flux_, args.x, kt2, form_fac_.get(), str_fun_.get(), args.mi2, args.mf2) / kt2;
          }));
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
        if (x == 0. || !x_valid_range.contains(x))
          return 0.;
        const FluxArguments params{x, mp2_, mx * mx};
        return 2. * M_PI * integr_->integrate(*func_, params, q2_range_) / x;
      }

    private:
      const Beam::KTFlux flux_;
      const HeavyIon hi_;
      std::unique_ptr<formfac::Parameterisation> form_fac_;
      std::unique_ptr<strfun::Parameterisation> str_fun_;
      std::unique_ptr<utils::Function1D> func_;
      std::unique_ptr<AnalyticIntegrator> integr_;

      struct FluxArguments {
        double x{0.}, mi2{0.}, mf2{0.};
      };
    };
  }  // namespace collflux
}  // namespace cepgen
typedef cepgen::collflux::GammaIntegrated CF_GI;
REGISTER_COLLFLUX("GammaIntegrated", CF_GI)
