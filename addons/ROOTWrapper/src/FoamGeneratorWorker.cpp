/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

using namespace cepgen;
using namespace std::string_literals;

/// Foam generator worker algorithm as developed by S. Jadach
/// (Institute of Nuclear Physics, Krakow, PL)
class FoamGeneratorWorker final : public GeneratorWorker, public TFoamIntegrand {
public:
  explicit FoamGeneratorWorker(const ParametersList& params) : GeneratorWorker(params) {
    if (const auto& random_number_mode = steer<std::string>("rngEngine"s); random_number_mode == "Ranlux"s)
      random_number_generator_.reset(new TRandom1);
    else if (random_number_mode == "generic"s)
      random_number_generator_.reset(new TRandom2);
    else if (random_number_mode == "MersenneTwister"s)
      random_number_generator_.reset(new TRandom3);
    else
      throw CG_FATAL("FoamGeneratorWorker") << "Unrecognised random generator: \"" << random_number_mode << "\".";
    random_number_generator_->SetSeed(steer<unsigned long long>("seed"));

    //--- a bit of printout for debugging
    CG_WARNING("FoamGeneratorWorker") << "This wrapping of the Foam generation algorithm implemented in ROOT "
                                         "libraries is still experimental! Please use with care...";
  }

  static ParametersDescription description() {
    auto desc = GeneratorWorker::description();
    desc.setDescription("Foam generator worker");
    desc.add("rngEngine"s, "MersenneTwister"s)
        .allow("Ranlux"s)
        .allow("generic"s)
        .allow("MersenneTwister"s)
        .setDescription("random number generator engine");
    desc.add("nCalls"s, 100'000).setDescription("number of calls for the cell evaluation");
    desc.add("nCells"s, 1000);
    desc.add("nSampl"s, 200);
    desc.add("nBin"s, 8);
    desc.add("EvPerBin"s, 25);
    desc.add("verbosity"s, 0).setDescription("Verbosity level");
    desc.add("seed"s, 42ull);
    return desc;
  }

  void initialise() override {
    foam_ = std::make_unique<TFoam>("Foam");
    foam_->SetPseRan(random_number_generator_.get());
    foam_->SetnCells(steer<int>("nCells"s));
    foam_->SetnSampl(steer<int>("nSampl"s));
    foam_->SetnBin(steer<int>("nBin"s));
    foam_->SetEvPerBin(steer<int>("EvPerBin"s));
    foam_->SetChat(std::max(steer<int>("verbosity"s), 0));
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
  double Density(int num_dimensions, double* coordinates) override {
    if (integrand_)
      return integrand_->eval(std::vector(coordinates, coordinates + num_dimensions));
    throw CG_FATAL("FoamGeneratorWorker:density") << "Integrand object was not initialised!";
  }

private:
  std::unique_ptr<TFoam> foam_;
  std::unique_ptr<TRandom> random_number_generator_;
};
REGISTER_GENERATOR_WORKER("Foam", FoamGeneratorWorker);
