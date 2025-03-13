/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/ProcessVariablesAnalyser.h"

using namespace cepgen;

/// Foam general-purpose integration algorithm
/// as developed by S. Jadach (Institute of Nuclear Physics, Krakow, PL)
class FoamIntegrator final : public Integrator, public TFoamIntegrand {
public:
  explicit FoamIntegrator(const ParametersList& params) : Integrator(params) {}

  static ParametersDescription description() {
    auto desc = Integrator::description();
    desc.setDescription("FOAM general purpose MC integrator");
    desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("root"));
    desc.add<int>("nCalls", 100'000).setDescription("number of calls for the cell evaluation");
    desc.add<int>("nCells", 1000);
    desc.add<int>("nSampl", 200);
    desc.add<int>("nBin", 8);
    desc.add<int>("EvPerBin", 25);
    desc.add<int>("verbose", 0).setDescription("Verbosity level");
    return desc;
  }

  Value integrate(Integrand& integrand) override {
    integrand_ = &integrand;
    std::unique_ptr<TFoam> foam(new TFoam("Foam"));
    CG_DEBUG("Integrator:integrate") << "FOAM integrator built\n\t"
                                     << "Version: " << foam->GetVersion() << ".";
    foam->SetPseRan(random_number_generator_->engine<TRandom>());
    foam->SetnCells(steer<int>("nCells"));
    foam->SetnSampl(steer<int>("nSampl"));
    foam->SetnBin(steer<int>("nBin"));
    foam->SetEvPerBin(steer<int>("EvPerBin"));
    foam->SetChat(std::max(verbosity_, 0));
    foam->SetRho(this);
    foam->SetkDim(integrand_->size());
    Integrator::checkLimits(*integrand_);
    coordinates_.resize(integrand_->size());
    foam->Initialize();
    std::unique_ptr<utils::ProcessVariablesAnalyser> analyser;
    if (integrand.hasProcess())
      analyser.reset(
          new utils::ProcessVariablesAnalyser(dynamic_cast<ProcessIntegrand&>(integrand).process(), ParametersList{}));
    const auto num_calls = steer<int>("nCalls");
    for (int i = 0; i < num_calls; ++i) {
      foam->MakeEvent();
      if (analyser)
        analyser->feed(foam->GetMCwt() / num_calls);
    }
    if (analyser)
      analyser->analyse();
    //--- launch integration
    double norm, err;
    foam->Finalize(norm, err);

    double result, abs_error;
    foam->GetIntegMC(result, abs_error);
    for (const auto& lim : limits_) {
      result *= lim.range();
      abs_error *= lim.range();
    }
    const auto res = Value{result, abs_error};

    CG_DEBUG("FoamIntegrator").log([&](auto& log) {
      double eps = 5.e-4, avewt, wtmax, sigma;
      foam->GetWtParams(eps, avewt, wtmax, sigma);
      const double ncalls = foam->GetnCalls();
      const double effic = wtmax > 0 ? avewt / wtmax : 0.;
      log << "Result: " << res << "\n\t"
          << "Relative error: " << res.relativeUncertainty() * 100. << "%\n\t"
          << "Dispersion/<wt> = " << sigma << ", <wt> = " << avewt << ", <wt>/wtmax = " << effic << ",\n\t"
          << " for epsilon = " << eps << "\n\t"
          << " nCalls (initialisation only)= " << ncalls << ".";
    });
    return res;
  }

  /// Compute the weight for a given phase space point
  inline double Density(int num_dimensions, double* x) override {
    if (!integrand_)
      throw CG_FATAL("FoamDensity") << "Integrand object not yet initialised!";
    for (int i = 0; i < num_dimensions; ++i)
      coordinates_[i] = limits_.at(i).x(x[i]);
    return integrand_->eval(coordinates_);
  }

private:
  Integrand* integrand_{nullptr};
  std::vector<double> coordinates_;
};
REGISTER_INTEGRATOR("Foam", FoamIntegrator);
