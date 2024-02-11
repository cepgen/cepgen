/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/TimeKeeper.h"
#include "CepGenAddOns/BasesWrapper/BasesCommonBlocks.h"

namespace cepgen {
  class SpringGeneratorWorker final : public GeneratorWorker {
  public:
    explicit SpringGeneratorWorker(const ParametersList& params)
        : GeneratorWorker(params), max_trials_(steer<int>("maxTrials")) {
      bscntl_.ipnt = steer<int>("verbose");
    }
    virtual ~SpringGeneratorWorker() {
      int lu = 6;
      spinfo_(lu);
    }

    void initialise() override {
      gIntegrand = dynamic_cast<Integrand*>(integrand_.get());
      sprng2_.mxtryp = max_trials_;
      sprng2_.nevent = 0;
    }
    bool next() override {
      if (!integrator_)
        throw CG_FATAL("SpringGeneratorWorker:next") << "No integrator object handled!";
      if (integrator_->name() != "bases")
        throw CG_FATAL("SpringGeneratorWorker:next") << "Spring generator is only compatible with Bases integrator.";

      CG_TICKER(const_cast<RunParameters*>(params_)->timeKeeper());

      sprng2_.ntrial = 0;
      sprng2_.miss = 0;
      int mxtry = max_trials_;
      spring_(integrand_call, mxtry);
      if (sprng2_.miss)
        return false;

      // return with an accepted event
      return storeEvent();
    }

    static ParametersDescription description() {
      auto desc = GeneratorWorker::description();
      desc.setDescription("Spring/Bases worker");
      desc.add<int>("maxTrials", 50).setDescription("maximum number of trials per generation");
      desc.add<int>("verbose", 0);
      return desc;
    }

  private:
    static Integrand* gIntegrand;
    static double integrand_call(double in[]) {
      if (!gIntegrand)
        throw CG_FATAL("SpringGeneratorWorker") << "Integrand was not specified before event generation.";
      return gIntegrand->eval(std::vector<double>(in, in + gIntegrand->size()));
    }

    const int max_trials_;
  };
  Integrand* SpringGeneratorWorker::gIntegrand = nullptr;
}  // namespace cepgen
REGISTER_GENERATOR_WORKER("spring", SpringGeneratorWorker);
