/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
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

#include <TFoam.h>
#include <TFoamIntegrand.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/ProcessVariablesAnalyser.h"

namespace cepgen {
  /// Foam general-purpose integration algorithm
  /// as developed by S. Jadach (Institute of Nuclear Physics, Krakow, PL)
  class IntegratorFoam final : public Integrator, public TFoamIntegrand {
  public:
    explicit IntegratorFoam(const ParametersList&);

    static ParametersDescription description();

    Value integrate(Integrand&) override;

    /// Compute the weight for a given phase space point
    inline double Density(int ndim, double* x) override {
      if (!integrand_)
        throw CG_FATAL("FoamDensity") << "Integrand object not yet initialised!";
      for (int i = 0; i < ndim; ++i)
        coord_[i] = limits_.at(i).x(x[i]);
      return integrand_->eval(coord_);
    }

  private:
    std::unique_ptr<TFoam> foam_;
    Integrand* integrand_{nullptr};
    std::vector<double> coord_;
  };

  IntegratorFoam::IntegratorFoam(const ParametersList& params) : Integrator(params), foam_(new TFoam("Foam")) {
    CG_DEBUG("Integrator:build") << "FOAM integrator built\n\t"
                                 << "Version: " << foam_->GetVersion() << ".";
  }

  Value IntegratorFoam::integrate(Integrand& integrand) {
    integrand_ = &integrand;
    foam_.reset(new TFoam("Foam"));
    foam_->SetPseRan(rnd_gen_->engine<TRandom>());
    foam_->SetnCells(steer<int>("nCells"));
    foam_->SetnSampl(steer<int>("nSampl"));
    foam_->SetnBin(steer<int>("nBin"));
    foam_->SetEvPerBin(steer<int>("EvPerBin"));
    foam_->SetChat(std::max(verbosity_, 0));
    foam_->SetRho(this);
    foam_->SetkDim(integrand_->size());
    Integrator::checkLimits(*integrand_);
    coord_.resize(integrand_->size());
    foam_->Initialize();
    std::unique_ptr<utils::ProcessVariablesAnalyser> analyser;
    if (integrand.hasProcess())
      analyser.reset(
          new utils::ProcessVariablesAnalyser(dynamic_cast<ProcessIntegrand&>(integrand).process(), ParametersList{}));
    const auto num_calls = steer<int>("nCalls");
    for (int i = 0; i < num_calls; ++i) {
      foam_->MakeEvent();
      if (analyser)
        analyser->feed(foam_->GetMCwt() / num_calls);
    }
    if (analyser)
      analyser->analyse();
    //--- launch integration
    double norm, err;
    foam_->Finalize(norm, err);

    double result, abs_error;
    foam_->GetIntegMC(result, abs_error);
    for (const auto& lim : limits_) {
      result *= lim.range();
      abs_error *= lim.range();
    }
    const auto res = Value{result, abs_error};

    CG_DEBUG("IntegratorFoam").log([&](auto& log) {
      double eps = 5.e-4, avewt, wtmax, sigma;
      foam_->GetWtParams(eps, avewt, wtmax, sigma);
      const double ncalls = foam_->GetnCalls();
      const double effic = wtmax > 0 ? avewt / wtmax : 0.;
      log << "Result: " << res << "\n\t"
          << "Relative error: " << res.relativeUncertainty() * 100. << "%\n\t"
          << "Dispersion/<wt> = " << sigma << ", <wt> = " << avewt << ", <wt>/wtmax = " << effic << ",\n\t"
          << " for epsilon = " << eps << "\n\t"
          << " nCalls (initialisation only)= " << ncalls << ".";
    });
    return res;
  }

  ParametersDescription IntegratorFoam::description() {
    auto desc = Integrator::description();
    desc.setDescription("FOAM general purpose MC integrator");
    desc.add<ParametersDescription>("randomGenerator", ParametersDescription().setName<std::string>("root"));
    desc.add<int>("nCalls", 100'000).setDescription("number of calls for the cell evaluation");
    desc.add<int>("nCells", 1000);
    desc.add<int>("nSampl", 200);
    desc.add<int>("nBin", 8);
    desc.add<int>("EvPerBin", 25);
    desc.add<int>("verbose", 0).setDescription("Verbosity level");
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("Foam", IntegratorFoam);
