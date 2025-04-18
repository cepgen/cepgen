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
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

/// Foam general-purpose integration algorithm
/// as developed by S. Jadach (Institute of Nuclear Physics, Krakow, PL)
class FoamIntegrator final : public Integrator, public TFoamIntegrand {
public:
  explicit FoamIntegrator(const ParametersList& params)
      : Integrator(params),
        random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {}

  static ParametersDescription description() {
    auto desc = Integrator::description();
    desc.setDescription("FOAM general purpose MC integrator");
    desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("root"));
    desc.add("nCalls", 100'000).setDescription("number of calls for the cell evaluation");
    desc.add("nCells", 1000).setDescription("number of allocated number of cells");
    desc.add("nSampl", 200).setDescription("number of MC events in the cell MC exploration");
    desc.add("nBin", 8).setDescription("number of bins in edge-histogram in cell exploration");
    desc.add("OptRej", 1)
        .allow(0, "weighted events")
        .allow(1, "unweighted events")
        .setDescription("MC events weight determination type");
    desc.add("OptDrive", 2).setDescription("maximum weight reduction (1 for variance reduction)");
    desc.add("MaxWtRej", 1.1).setDescription("maximum weight used to get unweighted MC events");
    desc.add("EvPerBin", 25)
        .setDescription("maximum number of the effective wt=1 events/bin (0 deactivates this option)");
    return desc;
  }

  Value run(Integrand& integrand, const std::vector<Limits>& range) override {
    integrand_ = &integrand;
    range_ = range;
    std::unique_ptr<TFoam> foam(new TFoam("Foam"));
    CG_DEBUG("Integrator:integrate") << "FOAM integrator built\n\t"
                                     << "Version: " << foam->GetVersion() << ".";
    foam->SetPseRan(random_generator_->engine<TRandom>());
    foam->SetnCells(steer<int>("nCells"));
    foam->SetnSampl(steer<int>("nSampl"));
    foam->SetnBin(steer<int>("nBin"));
    foam->SetOptRej(steer<int>("OptRej"));
    foam->SetOptDrive(steer<int>("OptDrive"));
    foam->SetMaxWtRej(steer<double>("MaxWtRej"));
    foam->SetEvPerBin(steer<int>("EvPerBin"));
    foam->SetChat(std::max(verbosity_, 0));
    foam->SetRho(this);
    foam->SetkDim(integrand_->size());
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
    for (const auto& lim : range) {
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
    std::vector<double> coordinates;
    for (int i = 0; i < num_dimensions; ++i)
      coordinates.emplace_back(range_.at(i).x(x[i]));
    return integrand_->eval(coordinates);
  }

private:
  const std::unique_ptr<utils::RandomGenerator> random_generator_;
  Integrand* integrand_{nullptr};
  std::vector<Limits> range_;
};
REGISTER_INTEGRATOR("Foam", FoamIntegrator);
