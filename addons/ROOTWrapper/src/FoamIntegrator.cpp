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
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/ProcessVariablesAnalyser.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace std::string_literals;

/// Foam integration algorithm
/// as developed by S. Jadach (Institute of Nuclear Physics, Krakow, PL)
class FoamIntegrator final : public Integrator, public TFoamIntegrand {
public:
  explicit FoamIntegrator(const ParametersList& params)
      : Integrator(params),
        random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"s))) {}

  static ParametersDescription description() {
    auto desc = Integrator::description();
    desc.setDescription("FOAM general purpose MC integrator");
    desc.add("randomGenerator"s, RandomGeneratorFactory::get().describeParameters("root"));
    desc.add("nCalls"s, 100'000).setDescription("number of calls for the cell evaluation");
    desc.add("nCells"s, 1000).setDescription("number of allocated number of cells");
    desc.add("nSampl"s, 200).setDescription("number of MC events in the cell MC exploration");
    desc.add("nBin"s, 8).setDescription("number of bins in edge-histogram in cell exploration");
    desc.add("OptRej"s, 1)
        .allow(0, "weighted events")
        .allow(1, "unweighted events")
        .setDescription("MC events weight determination type");
    desc.add("OptDrive"s, 2).setDescription("maximum weight reduction (1 for variance reduction)");
    desc.add("MaxWtRej"s, 1.1).setDescription("maximum weight used to get unweighted MC events");
    desc.add("EvPerBin"s, 25)
        .setDescription("maximum number of the effective wt=1 events/bin (0 deactivates this option)");
    return desc;
  }

  Value run(Integrand& integrand, const std::vector<Limits>& range) override {
    integrand_ = &integrand;
    range_ = range;
    const auto foam = std::make_unique<TFoam>("Foam");
    CG_DEBUG("Integrator:integrate") << "FOAM integrator built\n\t"
                                     << "Version: " << foam->GetVersion() << ".";
    foam->SetPseRan(random_generator_->engine<TRandom>());
    foam->SetnCells(steer<int>("nCells"s));
    foam->SetnSampl(steer<int>("nSampl"s));
    foam->SetnBin(steer<int>("nBin"s));
    foam->SetOptRej(steer<int>("OptRej"s));
    foam->SetOptDrive(steer<int>("OptDrive"s));
    foam->SetMaxWtRej(steer<double>("MaxWtRej"s));
    foam->SetEvPerBin(steer<int>("EvPerBin"s));
    foam->SetChat(std::max(verbosity_, 0));
    foam->SetRho(this);
    foam->SetkDim(integrand_->size());
    foam->Initialize();
    std::unique_ptr<utils::ProcessVariablesAnalyser> analyser;
    if (integrand.hasProcess())
      analyser.reset(
          new utils::ProcessVariablesAnalyser(dynamic_cast<ProcessIntegrand&>(integrand).process(), ParametersList{}));
    const auto num_calls = steer<int>("nCalls");
    // launch integration
    for (int i = 0; i < num_calls; ++i) {
      foam->MakeEvent();
      if (analyser)
        analyser->feed(foam->GetMCwt() / num_calls);
    }
    if (analyser)
      analyser->analyse();
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
      double eps = 5.e-4, average_event_weight, maximum_weight, sigma;
      foam->GetWtParams(eps, average_event_weight, maximum_weight, sigma);
      const double num_function_calls = foam->GetnCalls();
      const double efficiency = maximum_weight > 0 ? average_event_weight / maximum_weight : 0.;
      log << "Result: " << res << "\n\t"
          << "Relative error: " << res.relativeUncertainty() * 100. << "%\n\t"
          << "Dispersion/<wt> = " << sigma << ", <wt> = " << average_event_weight << ", <wt>/wtmax = " << efficiency
          << ",\n\t"
          << " for epsilon = " << eps << "\n\t"
          << " nCalls (initialisation only)= " << num_function_calls << ".";
    });
    return res;
  }

  /// Compute the weight for a given phase space point
  double Density(int num_dimensions, double* coordinates) override {
    if (!integrand_)
      throw CG_FATAL("FoamDensity") << "Integrand object not yet initialised!";
    std::vector<double> vec_coordinates;
    for (int i = 0; i < num_dimensions; ++i)
      vec_coordinates.emplace_back(range_.at(i).x(coordinates[i]));
    return integrand_->eval(vec_coordinates);
  }

private:
  const std::unique_ptr<utils::RandomGenerator> random_generator_;
  Integrand* integrand_{nullptr};
  std::vector<Limits> range_;
};
REGISTER_INTEGRATOR("Foam", FoamIntegrator);
