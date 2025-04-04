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
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

using namespace cepgen;

class GridOptimisedGeneratorWorker final : public GeneratorWorker {
public:
  /// Book the memory slots and structures for the generator
  explicit GridOptimisedGeneratorWorker(const ParametersList& params) : GeneratorWorker(params) {}

  void initialise() override;
  bool next() override;

  static ParametersDescription description() {
    auto desc = GeneratorWorker::description();
    desc.setDescription("Grid-optimised worker");
    desc.add("binSize", 3);
    return desc;
  }

private:
  static constexpr int UNASSIGNED_BIN = -999;  ///< Placeholder for invalid bin indexing

  bool correctionCycle(bool&);               ///< Apply a correction cycle to the grid
  void computeGenerationParameters() const;  ///< Prepare the object for event generation

  std::unique_ptr<GridParameters> grid_;  ///< Set of parameters for the integration/event generation grid
  int ps_bin_{UNASSIGNED_BIN};            ///< Last bin to be corrected
  std::vector<double> coords_;            ///< Phase space coordinates being evaluated
};

void GridOptimisedGeneratorWorker::initialise() {
  grid_ = std::make_unique<GridParameters>(steer<int>("binSize"), integrand_->size());
  coords_ = std::vector<double>(integrand_->size());
  if (!grid_->prepared())
    computeGenerationParameters();
  CG_DEBUG("GridOptimisedGeneratorWorker:initialise")
      << "Dim-" << integrand_->size() << " " << integrator_->name() << " integrator "
      << "set for dim-" << grid_->n(0).size() << " grid.";
}

//-----------------------------------------------------------------------------------------------
// events generation part
//-----------------------------------------------------------------------------------------------

bool GridOptimisedGeneratorWorker::next() {
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
      ps_bin_ = integrator_->uniform({0., static_cast<double>(grid_->size())});
      y = integrator_->uniform({0., grid_->globalMax()});
      grid_->increment(ps_bin_);
    } while (y > grid_->maxValue(ps_bin_));
    // shoot a point x in this bin
    grid_->shoot(integrator_, ps_bin_, coords_);
    if (weight = integrator_->eval(*integrand_, coords_);  // get weight for selected x value
        weight > y)
      break;
  }

  if (weight > grid_->maxValue(ps_bin_))
    // if weight is higher than local or global maximum, init correction cycle for the next event
    grid_->initCorrectionCycle(ps_bin_, weight);
  else  // no grid correction needed for this bin
    ps_bin_ = UNASSIGNED_BIN;

  // return with an accepted event
  return storeEvent();
}

bool GridOptimisedGeneratorWorker::correctionCycle(bool& store) {
  CG_TICKER(const_cast<RunParameters*>(run_params_)->timeKeeper());

  CG_DEBUG_LOOP("GridOptimisedGeneratorWorker:correction") << "Correction cycles are started.\n\t"
                                                           << "bin = " << ps_bin_ << "\n\t"
                                                           << "correction value = " << grid_->correctionValue() << ".";

  if (grid_->correctionValue() >= 1.)
    grid_->setCorrectionValue(grid_->correctionValue() - 1.);

  if (integrator_->uniform() < grid_->correctionValue()) {
    grid_->setCorrectionValue(-1.);
    grid_->shoot(integrator_, ps_bin_, coords_);  // select x values in phase space bin
    const double weight = integrator_->eval(*integrand_, coords_);
    grid_->rescale(ps_bin_, weight);  // parameter for correction of correction
    if (weight >= integrator_->uniform({0., grid_->maxValueDiff()}) + grid_->maxHistValue()) {
      store = true;  // accept event
      return true;
    }
    return false;
  }
  // correction if too big weight is found while correction
  // (all your bases are belong to us...)
  return grid_->correct(ps_bin_);
}

//-----------------------------------------------------------------------------------------------
// initial preparation run before the generation of unweighted events
//-----------------------------------------------------------------------------------------------

void GridOptimisedGeneratorWorker::computeGenerationParameters() const {
  if (!run_params_)
    throw CG_FATAL("GridOptimisedGeneratorWorker:setGen") << "No steering parameters specified!";
  if (!integrator_)
    throw CG_FATAL("GridOptimisedGeneratorWorker:setGen") << "No integrator object specified!";

  integrand_->setStorage(false);

  CG_INFO("GridOptimisedGeneratorWorker:setGen")
      << "Preparing the grid (" << utils::s("point", run_params_->generation().numPoints(), true) << "/bin) "
      << "for the generation of unweighted events.";

  const double inv_num_points = 1. / run_params_->generation().numPoints();
  std::vector point_coord(integrand_->size(), 0.);
  if (point_coord.size() < grid_->n(0).size())
    throw CG_FATAL("GridParameters:shoot") << "Coordinates vector multiplicity is insufficient!";

  // ...
  double sum = 0., sum2 = 0., sum2p = 0.;

  utils::ProgressBar progress_bar(grid_->size(), 5);

  // main loop
  for (unsigned int i = 0; i < grid_->size(); ++i) {
    double fsum = 0., fsum2 = 0.;
    for (size_t j = 0; j < run_params_->generation().numPoints(); ++j) {
      grid_->shoot(integrator_, i, point_coord);
      const double weight = integrator_->eval(*integrand_, point_coord);
      grid_->setValue(i, weight);
      fsum += weight;
      fsum2 += weight * weight;
    }
    const double av = fsum * inv_num_points, av2 = fsum2 * inv_num_points;
    const double sig2 = av2 - av * av;
    sum += av;
    sum2 += av2;
    sum2p += sig2;

    // per-bin debugging loop
    CG_DEBUG_LOOP("GridOptimisedGeneratorWorker:setGen").log([&](auto& dbg) {
      const double sig = sqrt(sig2);
      const double eff = (grid_->maxValue(i) != 0.) ? av / grid_->maxValue(i) : 0.;
      dbg << "n-vector for bin " << i << ": " << utils::repr(grid_->n(i)) << "\n\t"
          << "av   = " << av << "\n\t"
          << "sig  = " << sig << "\n\t"
          << "fmax = " << grid_->maxValue(i) << "\n\t"
          << "eff  = " << eff;
    });
    progress_bar.update(i + 1);
  }  // end of main loop

  const double inv_max = 1. / grid_->size();
  sum *= inv_max;
  sum2 *= inv_max;
  sum2p *= inv_max;

  const double sig = sqrt(sum2 - sum * sum), sigp = sqrt(sum2p);

  double eff1 = 0.;
  for (unsigned int i = 0; i < grid_->size(); ++i)
    eff1 += sum / grid_->size() * grid_->maxValue(i);
  const double eff2 = sum / grid_->globalMax();

  CG_DEBUG("GridOptimisedGeneratorWorker:setGen") << "Average function value         = " << sum << "\n\t"
                                                  << "Average squared function value = " << sum2 << "\n\t"
                                                  << "Overall standard deviation     = " << sig << "\n\t"
                                                  << "Average standard deviation     = " << sigp << "\n\t"
                                                  << "Maximum function value         = " << grid_->globalMax() << "\n\t"
                                                  << "Average inefficiency           = " << eff1 << "\n\t"
                                                  << "Overall inefficiency           = " << eff2;
  grid_->setPrepared(true);
  integrand_->setStorage(true);  // from now on events will be stored

  CG_INFO("GridOptimisedGeneratorWorker:setGen")
      << "Finished the grid preparation. Now launching the unweighted event production.";
}
REGISTER_GENERATOR_WORKER("grid_optimised", GridOptimisedGeneratorWorker);
