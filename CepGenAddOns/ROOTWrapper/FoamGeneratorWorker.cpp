/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/GeneratorWorker.h"
#include "CepGen/Integration/ProcessIntegrand.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"

namespace cepgen {
  /// Foam generator worker algorithm as developed by S. Jadach
  /// (Institute of Nuclear Physics, Krakow, PL)
  class FoamGeneratorWorker final : public GeneratorWorker, public TFoamIntegrand {
  public:
    explicit FoamGeneratorWorker(const ParametersList& params) : GeneratorWorker(params) {
      const auto& rnd_mode = steer<std::string>("rngEngine");
      if (rnd_mode == "Ranlux")
        rnd_.reset(new TRandom1);
      else if (rnd_mode == "generic")
        rnd_.reset(new TRandom2);
      else if (rnd_mode == "MersenneTwister")
        rnd_.reset(new TRandom3);
      else
        throw CG_FATAL("FoamGeneratorWorker") << "Unrecognised random generator: \"" << rnd_mode << "\".";
      rnd_->SetSeed(steer<unsigned long long>("seed"));

      //--- a bit of printout for debugging
      CG_WARNING("FoamGeneratorWorker") << "This wrapping of the Foam generation algorithm implemented in ROOT "
                                           "libraries is still experimental! Please use with care...";
    }

    static ParametersDescription description() {
      auto desc = GeneratorWorker::description();
      desc.setDescription("Foam generator worker");
      desc.add<std::string>("rngEngine", "MersenneTwister")
          .setDescription(
              "Set random number generator engine (currently handled: 'Ranlux', 'generic', 'MersenneTwister')");
      desc.add<int>("nCalls", 100'000).setDescription("number of calls for the cell evaluation");
      desc.add<int>("nCells", 1000);
      desc.add<int>("nSampl", 200);
      desc.add<int>("nBin", 8);
      desc.add<int>("EvPerBin", 25);
      desc.add<int>("verbose", 0).setDescription("Verbosity level");
      desc.add<unsigned long long>("seed", 42ull);
      return desc;
    }

    std::unique_ptr<GeneratorWorker> clone() const override {
      return std::unique_ptr<GeneratorWorker>(new FoamGeneratorWorker(parameters()));
    }

    void initialise() override {
      foam_.reset(new TFoam("Foam"));
      foam_->SetPseRan(rnd_.get());
      foam_->SetnCells(steer<int>("nCells"));
      foam_->SetnSampl(steer<int>("nSampl"));
      foam_->SetnBin(steer<int>("nBin"));
      foam_->SetEvPerBin(steer<int>("EvPerBin"));
      foam_->SetChat(std::max(steer<int>("verbose"), 0));
      foam_->SetRho(this);
      foam_->SetkDim(integrand_->size());
      foam_->Initialize();
      CG_DEBUG("FoamGeneratorWorker:build") << "FOAM integrator built\n\t"
                                            << "Version: " << foam_->GetVersion() << ".";
    }
    bool next() override {
      foam_->MakeEvent();
      return storeEvent();
    }

    /// Compute the weight for a given phase space point
    inline double Density(int ndim, double* x) override {
      if (integrand_)
        return integrand_->eval(std::vector<double>(x, x + ndim));
      throw CG_FATAL("FoamGeneratorWorker:density") << "Integrand object was not initialised!";
    }

  private:
    std::unique_ptr<TFoam> foam_;
    std::unique_ptr<TRandom> rnd_;
    std::vector<double> coord_;
  };
}  // namespace cepgen

REGISTER_GENERATOR_WORKER("Foam", FoamGeneratorWorker);
