/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Integration/GridParameters.h"
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ProgressBar.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  class TrivialGeneratorWorker final : public GeneratorWorker {
  public:
    /// Book the memory slots and structures for the generator
    explicit TrivialGeneratorWorker(const ParametersList& params) : GeneratorWorker(params) {}

    void initialise() override;
    bool next() override;

    static ParametersDescription description() {
      auto desc = GeneratorWorker::description();
      desc.setDescription("Grid-optimised worker");
      desc.add<int>("binSize", 3);
      return desc;
    }

  private:
    /// Placeholder for invalid bin indexing
    static constexpr int UNASSIGNED_BIN = -999;

    /// Apply a correction cycle to the grid
    bool correctionCycle(bool&);
    /// Prepare the object for event generation
    void computeGenerationParameters();

    /// Set of parameters for the integration/event generation grid
    std::unique_ptr<GridParameters> grid_;
    /// Selected bin at which the function will be evaluated
    int ps_bin_{UNASSIGNED_BIN};  ///< Last bin to be corrected
    std::vector<double> coords_;  ///< Phase space coordinates being evaluated
  };

  void TrivialGeneratorWorker::initialise() {
    grid_.reset(new GridParameters(steer<int>("binSize"), integrand_->size()));
    coords_ = std::vector<double>(integrand_->size());
    if (!grid_->prepared())
      computeGenerationParameters();
    CG_DEBUG("TrivialGeneratorWorker:initialise")
        << "Dim-" << integrand_->size() << " " << integrator_->name() << " integrator "
        << "set for dim-" << grid_->n(0).size() << " grid.";
  }

  //-----------------------------------------------------------------------------------------------
  // events generation part
  //-----------------------------------------------------------------------------------------------

  bool TrivialGeneratorWorker::next() {
    if (!integrator_)
      throw CG_FATAL("TrivialGeneratorWorker:next") << "No integrator object handled!";
    if (!grid_)
      throw CG_FATAL("TrivialGeneratorWorker:next") << "Grid object was not initialised.";

    CG_TICKER(const_cast<Parameters*>(params_)->timeKeeper());

    // apply correction cycles if required from previous event
    if (ps_bin_ != UNASSIGNED_BIN) {
      bool store = false;
      while (!correctionCycle(store)) {
      }
      if (store)
        return storeEvent();
    }

    //--- normal generation cycle

    double weight = 0.;
    while (true) {
      double y = -1.;
      // select a function value and reject if fmax is too small
      do {
        // ...
        ps_bin_ = integrator_->uniform(0., grid_->size());
        y = integrator_->uniform(0., grid_->globalMax());
        grid_->increment(ps_bin_);
      } while (y > grid_->maxValue(ps_bin_));
      // shoot a point x in this bin
      grid_->shoot(integrator_, ps_bin_, coords_);
      // get weight for selected x value
      weight = integrator_->eval(*integrand_, coords_);
      if (weight > y)
        break;
    }

    if (weight > grid_->maxValue(ps_bin_))
      // if weight is higher than local or global maximum,
      // init correction cycle for the next event
      grid_->initCorrectionCycle(ps_bin_, weight);
    else  // no grid correction needed for this bin
      ps_bin_ = UNASSIGNED_BIN;

    // return with an accepted event
    return storeEvent();
  }

  bool TrivialGeneratorWorker::correctionCycle(bool& store) {
    CG_TICKER(const_cast<Parameters*>(params_)->timeKeeper());

    CG_DEBUG_LOOP("TrivialGeneratorWorker:correction") << "Correction cycles are started.\n\t"
                                                       << "bin = " << ps_bin_ << "\n\t"
                                                       << "correction value = " << grid_->correctionValue() << ".";

    if (grid_->correctionValue() >= 1.)
      grid_->setCorrectionValue(grid_->correctionValue() - 1.);

    if (integrator_->uniform() < grid_->correctionValue()) {
      grid_->setCorrectionValue(-1.);
      // select x values in phase space bin
      grid_->shoot(integrator_, ps_bin_, coords_);
      const double weight = integrator_->eval(*integrand_, coords_);
      // parameter for correction of correction
      grid_->rescale(ps_bin_, weight);
      // accept event
      if (weight >= integrator_->uniform(0., grid_->maxValueDiff()) + grid_->maxHistValue()) {
        store = true;
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

  void TrivialGeneratorWorker::computeGenerationParameters() {
    if (!params_)
      throw CG_FATAL("TrivialGeneratorWorker:setGen") << "No steering parameters specified!";
    if (!integrator_)
      throw CG_FATAL("TrivialGeneratorWorker:setGen") << "No integrator object specified!";

    integrand_->setStorage(false);

    CG_INFO("TrivialGeneratorWorker:setGen")
        << "Preparing the grid (" << utils::s("point", params_->generation().numPoints(), true) << "/bin) "
        << "for the generation of unweighted events.";

    const double inv_num_points = 1. / params_->generation().numPoints();
    std::vector<double> point_coord(integrand_->size(), 0.);
    if (point_coord.size() < grid_->n(0).size())
      throw CG_FATAL("GridParameters:shoot") << "Coordinates vector multiplicity is insufficient!";

    // ...
    double sum = 0., sum2 = 0., sum2p = 0.;

    utils::ProgressBar prog_bar(grid_->size(), 5);

    //--- main loop
    for (unsigned int i = 0; i < grid_->size(); ++i) {
      double fsum = 0., fsum2 = 0.;
      for (size_t j = 0; j < params_->generation().numPoints(); ++j) {
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
      CG_DEBUG_LOOP("TrivialGeneratorWorker:setGen").log([&](auto& dbg) {
        const double sig = sqrt(sig2);
        const double eff = (grid_->maxValue(i) != 0.) ? av / grid_->maxValue(i) : 0.;
        dbg << "n-vector for bin " << i << ": " << utils::repr(grid_->n(i)) << "\n\t"
            << "av   = " << av << "\n\t"
            << "sig  = " << sig << "\n\t"
            << "fmax = " << grid_->maxValue(i) << "\n\t"
            << "eff  = " << eff;
      });
      prog_bar.update(i + 1);
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

    CG_DEBUG("TrivialGeneratorWorker:setGen") << "Average function value         = " << sum << "\n\t"
                                              << "Average squared function value = " << sum2 << "\n\t"
                                              << "Overall standard deviation     = " << sig << "\n\t"
                                              << "Average standard deviation     = " << sigp << "\n\t"
                                              << "Maximum function value         = " << grid_->globalMax() << "\n\t"
                                              << "Average inefficiency           = " << eff1 << "\n\t"
                                              << "Overall inefficiency           = " << eff2;
    grid_->setPrepared(true);
    //--- from now on events will be stored
    integrand_->setStorage(true);

    CG_INFO("TrivialGeneratorWorker:setGen") << "Grid prepared! Now launching the production.";
  }
}  // namespace cepgen
REGISTER_GENERATOR_WORKER("trivial", TrivialGeneratorWorker);
