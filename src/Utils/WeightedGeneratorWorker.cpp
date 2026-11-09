/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

/// A worker producing weighted events, without any accept/reject step
///
/// Phase space points are drawn uniformly in the unit hypercube, and mapped by the integrator (e.g. through the
/// Vegas grid if its "treat" option is enabled). Event weights are normalised to a unit mean over non-zero events.
class WeightedGeneratorWorker final : public GeneratorWorker {
public:
  explicit WeightedGeneratorWorker(const ParametersList& params)
      : GeneratorWorker(params),
        random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {}

  static ParametersDescription description() {
    auto desc = GeneratorWorker::description();
    desc.setDescription("Weighted events worker");
    desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
        .setDescription("random number generator engine");
    desc.add("numPresamplingPoints", 100'000).setDescription("number of points used to normalise the event weights");
    return desc;
  }

  void initialise() override {
    if (!integrator_)
      throw CG_FATAL("WeightedGeneratorWorker:initialise") << "No integrator object handled!";
    coordinates_.resize(integrand_->size());

    integrand_->setStorage(false);  // no event storage while normalising the weights
    const auto num_points = steer<int>("numPresamplingPoints");
    auto sum_weights = 0.;
    size_t num_nonzero = 0;
    for (int i = 0; i < num_points; ++i)
      if (const auto weight = shoot(); weight > 0.)
        sum_weights += weight, ++num_nonzero;
    if (num_nonzero == 0)
      throw CG_FATAL("WeightedGeneratorWorker:initialise")
          << "None of the " << num_points << " presampling points passed the cuts.";
    mean_weight_ = sum_weights / num_nonzero;
    CG_INFO("WeightedGeneratorWorker:initialise")
        << "Presampling finished: " << num_nonzero << "/" << num_points << " non-zero points, "
        << "mean non-zero weight: " << mean_weight_ << " pb.";
    integrand_->setStorage(true);  // from now on events will be stored
  }

  bool next() override {
    double weight;
    do  // points failing the cuts have a null weight, and are never stored
      weight = shoot();
    while (weight <= 0.);
    return storeEvent(weight / mean_weight_);
  }
  bool weighted() const final { return true; }

private:
  /// Evaluate the integrand at a random point of the unit hypercube
  double shoot() {
    for (auto& coordinate : coordinates_)
      coordinate = random_generator_->uniform();
    return integrator_->eval(*integrand_, coordinates_);
  }

  const std::unique_ptr<utils::RandomGenerator> random_generator_;  ///< Random number generator
  std::vector<double> coordinates_;                                 ///< Phase space coordinates being evaluated
  double mean_weight_{1.};                                          ///< Mean non-zero integrand value, in pb
};
REGISTER_GENERATOR_WORKER("weighted", WeightedGeneratorWorker);
