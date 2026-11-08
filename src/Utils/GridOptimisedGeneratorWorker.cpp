/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen;

/// A Vegas grid-aware optimised event generator
class GridOptimisedGeneratorWorker final : public GeneratorWorker {
public:
  explicit GridOptimisedGeneratorWorker(const ParametersList& params)
      : GeneratorWorker(params),
        random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {}

  static ParametersDescription description() {
    auto desc = GeneratorWorker::description();
    desc.setDescription("Grid-optimised worker");
    desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
        .setDescription("random number generator engine");
    desc.add("binSize", 3);
    return desc;
  }

  void initialise() override {
    grid_ = std::make_unique<GridParameters>(steer<int>("binSize"), integrand_->size());
    coordinates_ = std::vector<double>(integrand_->size());
    if (!grid_->prepared())
      computeGenerationParameters();
    CG_DEBUG("GridOptimisedGeneratorWorker:initialise")
        << "Dim-" << integrand_->size() << " " << integrator_->name() << " integrator "
        << "set for dim-" << grid_->n(0).size() << " grid.";
  }
  bool next() override {
    if (!integrator_)
      throw CG_FATAL("GridOptimisedGeneratorWorker:next") << "No integrator object handled!";
    if (!grid_)
      throw CG_FATAL("GridOptimisedGeneratorWorker:next") << "Grid object was not initialised.";

    CG_TICKER(const_cast<RunParameters*>(run_params_)->timeKeeper());
    if (ps_bin_ != UNASSIGNED_BIN) {  // apply correction cycles if required from previous event
      bool store = false;
      while (!correctionCycle(store)) {
      }
      if (store)
        return storeEvent();
    }
    // normal generation cycle
    double weight;
    while (true) {
      double y;
      do {  // select a function value and reject if fmax is too small
        ps_bin_ = random_generator_->uniformInt(0, grid_->size() - 1);
        y = random_generator_->uniform(0., grid_->globalMax());
        grid_->increment(ps_bin_);
      } while (y > grid_->maxValue(ps_bin_));
      grid_->shoot(*random_generator_, ps_bin_, coordinates_);    // shoot a point x in this bin
      if (weight = integrator_->eval(*integrand_, coordinates_);  // get weight for selected x value
          weight > y)
        break;
    }
    if (weight > grid_->maxValue(ps_bin_))          // if weight is higher than local or global maximum,
      grid_->initCorrectionCycle(ps_bin_, weight);  // init correction cycle for the next event
    else                                            // no grid correction needed for this bin
      ps_bin_ = UNASSIGNED_BIN;
    return storeEvent();  // return with an accepted event
  }

private:
  static constexpr int UNASSIGNED_BIN = -999;  ///< Placeholder for invalid bin indexing

  /// Apply a correction cycle to the grid
  bool correctionCycle(bool& store) {
    CG_TICKER(const_cast<RunParameters*>(run_params_)->timeKeeper());

    CG_DEBUG_LOOP("GridOptimisedGeneratorWorker:correction")
        << "Correction cycles are started.\n\t"
        << "bin = " << ps_bin_ << "\n\t"
        << "correction value = " << grid_->correctionValue() << ".";

    if (grid_->correctionValue() >= 1.)
      grid_->setCorrectionValue(grid_->correctionValue() - 1.);

    if (random_generator_->uniform() < grid_->correctionValue()) {
      grid_->setCorrectionValue(-1.);
      grid_->shoot(*random_generator_, ps_bin_, coordinates_);  // select x values in phase space bin
      const auto weight = integrator_->eval(*integrand_, coordinates_);
      grid_->rescale(ps_bin_, weight);  // parameter for correction of correction
      if (weight >= random_generator_->uniform(0., grid_->maxValueDiff()) + grid_->maxHistValue()) {
        store = true;  // accept event
        return true;
      }
      return false;
    }
    // correction if too big weight is found while correction
    // (all your bases are belong to us...)
    return grid_->correct(ps_bin_);
  }
  /// Prepare the object for event generation
  void computeGenerationParameters() const {
    if (!run_params_)
      throw CG_FATAL("GridOptimisedGeneratorWorker:setGen") << "No steering parameters specified!";
    if (!integrator_)
      throw CG_FATAL("GridOptimisedGeneratorWorker:setGen") << "No integrator object specified!";

    integrand_->setStorage(false);  // first deactivate the storage of event while we are fiddling around with the grid
    CG_INFO("GridOptimisedGeneratorWorker:setGen")
        << "Preparing the grid (" << utils::s("point", run_params_->generation().numPoints(), true) << "/bin) "
        << "for the generation of unweighted events.";

    const auto inv_num_points = 1. / run_params_->generation().numPoints();
    std::vector point_coord(integrand_->size(), 0.);
    if (point_coord.size() < grid_->n(0).size())
      throw CG_FATAL("GridParameters:setGen") << "Coordinates vector multiplicity is insufficient!";

    // main preparation loop
    utils::ProgressBar progress_bar(grid_->size(), 5);
    auto sum = 0., sum2 = 0., sum2p = 0.;
    for (size_t i = 0; i < grid_->size(); ++i) {
      auto fsum = 0., fsum2 = 0.;
      for (size_t j = 0; j < run_params_->generation().numPoints(); ++j) {
        grid_->shoot(*random_generator_, i, point_coord);
        const auto weight = integrator_->eval(*integrand_, point_coord);
        grid_->setValue(i, weight);
        fsum += weight;
        fsum2 += weight * weight;
      }
      const auto av = fsum * inv_num_points, av2 = fsum2 * inv_num_points, sig2 = av2 - av * av;
      sum += av;
      sum2 += av2;
      sum2p += sig2;
      CG_DEBUG_LOOP("GridOptimisedGeneratorWorker:setGen")
          .log([this, &i, &av, &sig2](auto& log) {  // per-bin debugging loop
            const double sig = sqrt(sig2);
            const double eff = (grid_->maxValue(i) != 0.) ? av / grid_->maxValue(i) : 0.;
            log << "n-vector for bin " << i << ": " << utils::repr(grid_->n(i)) << "\n\t"
                << "av   = " << av << "\n\t"
                << "sig  = " << sig << "\n\t"
                << "fmax = " << grid_->maxValue(i) << "\n\t"
                << "eff  = " << eff;
          });
      progress_bar.update(i + 1);
    }  // end of main preparation loop

    CG_DEBUG("GridOptimisedGeneratorWorker:setGen").log([this, &sum, &sum2, &sum2p](auto& log) {
      const double inv_max = 1. / grid_->size();
      sum *= inv_max;
      sum2 *= inv_max;
      sum2p *= inv_max;
      const auto sig = std::sqrt(sum2 - sum * sum), sigp = std::sqrt(sum2p);
      double eff1 = 0.;
      for (unsigned int i = 0; i < grid_->size(); ++i)
        eff1 += sum / grid_->size() * grid_->maxValue(i);
      const double eff2 = sum / grid_->globalMax();
      log << "Average function value         = " << sum << "\n\t"
          << "Average squared function value = " << sum2 << "\n\t"
          << "Overall standard deviation     = " << sig << "\n\t"
          << "Average standard deviation     = " << sigp << "\n\t"
          << "Maximum function value         = " << grid_->globalMax() << "\n\t"
          << "Average inefficiency           = " << eff1 << "\n\t"
          << "Overall inefficiency           = " << eff2;
    });
    grid_->setPrepared(true);
    integrand_->setStorage(true);  // from now on events will be stored
    CG_INFO("GridOptimisedGeneratorWorker:setGen")
        << "Finished the grid preparation. Now launching the unweighted event production.";
  }

  const std::unique_ptr<utils::RandomGenerator> random_generator_;  ///< Random number generator for grid population
  std::unique_ptr<GridParameters> grid_;  ///< Set of parameters for the integration/event generation grid
  int ps_bin_{UNASSIGNED_BIN};            ///< Last bin to be corrected
  std::vector<double> coordinates_;       ///< Phase space coordinates being evaluated
};
REGISTER_GENERATOR_WORKER("grid_optimised", GridOptimisedGeneratorWorker);
