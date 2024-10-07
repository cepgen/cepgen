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

// engines
#include <CLHEP/Random/DRand48Engine.h>
#include <CLHEP/Random/DualRand.h>
#include <CLHEP/Random/Hurd160Engine.h>
#include <CLHEP/Random/Hurd288Engine.h>
#include <CLHEP/Random/JamesRandom.h>
#include <CLHEP/Random/MTwistEngine.h>
#include <CLHEP/Random/MixMaxRng.h>
#include <CLHEP/Random/NonRandomEngine.h>
#include <CLHEP/Random/RandEngine.h>
#include <CLHEP/Random/RanecuEngine.h>
#include <CLHEP/Random/Ranlux64Engine.h>
#include <CLHEP/Random/RanluxEngine.h>
#include <CLHEP/Random/RanshiEngine.h>
#include <CLHEP/Random/TripleRand.h>
// distributions
#include <CLHEP/Random/RandBreitWigner.h>
#include <CLHEP/Random/RandExponential.h>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandGauss.h>
#include <CLHEP/Random/RandLandau.h>
#include <CLHEP/Random/RandPoisson.h>

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen::clhep {
  class RandomGenerator : public utils::RandomGenerator {
  public:
    explicit RandomGenerator(const ParametersList& params) : utils::RandomGenerator(params) {
      const auto& type = steer<std::string>("type");
      if (type == "HepJamesRandom")
        engine_ = std::make_unique<CLHEP::HepJamesRandom>();
      else if (type == "RandEngine")
        engine_ = std::make_unique<CLHEP::RandEngine>();
      else if (type == "DRand48Engine")
        engine_ = std::make_unique<CLHEP::DRand48Engine>();
      else if (type == "RanluxEngine")
        engine_ = std::make_unique<CLHEP::RanluxEngine>();
      else if (type == "Ranlux64Engine")
        engine_ = std::make_unique<CLHEP::Ranlux64Engine>();
      else if (type == "RanecuEngine")
        engine_ = std::make_unique<CLHEP::RanecuEngine>();
      else if (type == "Hurd160Engine")
        engine_ = std::make_unique<CLHEP::Hurd160Engine>();
      else if (type == "Hurd288Engine")
        engine_ = std::make_unique<CLHEP::Hurd288Engine>();
      else if (type == "MTwistEngine")
        engine_ = std::make_unique<CLHEP::MTwistEngine>();
      else if (type == "RanshiEngine")
        engine_ = std::make_unique<CLHEP::RanshiEngine>();
      else if (type == "DualRand")
        engine_ = std::make_unique<CLHEP::DualRand>();
      else if (type == "TripleRand")
        engine_ = std::make_unique<CLHEP::TripleRand>();
      else if (type == "NonRandomEngine")
        engine_ = std::make_unique<CLHEP::NonRandomEngine>();
      else
        throw CG_FATAL("clhep:RandomGenerator") << "Random number generator engine invalid: '" << type << "'.";

      engine_->setSeed(seed_, 0);
    }

    static ParametersDescription description() {
      auto desc = utils::RandomGenerator::description();
      desc.setDescription("CLHEP random number generator engine");
      desc.add<std::string>("type", "HepJamesRandom")
          .setDescription("random number engine")
          .allow("HepJamesRandom")
          .allow("RandEngine")
          .allow("DRand48Engine")
          .allow("RanluxEngine")
          .allow("Ranlux64Engine")
          .allow("RanecuEngine")
          .allow("Hurd160Engine")
          .allow("Hurd288Engine")
          .allow("MTwistEngine")
          .allow("RanshiEngine")
          .allow("DualRand")
          .allow("TripleRand")
          .allow("NonRandomEngine");
      return desc;
    }

    int uniformInt(int min, int max) override { return CLHEP::RandFlat::shootInt(engine_.get(), min, max + 1); }
    double uniform(double min, double max) override { return CLHEP::RandFlat::shoot(engine_.get(), min, max); }
    double normal(double mean, double rms) override { return CLHEP::RandGauss::shoot(engine_.get(), mean, rms); }
    double exponential(double exponent) override { return CLHEP::RandExponential::shoot(engine_.get(), exponent); }
    double breitWigner(double mean, double scale) override {
      return CLHEP::RandBreitWigner::shoot(engine_.get(), mean, scale);
    }
    double landau(double location, double width) override {
      return location + width * CLHEP::RandLandau::shoot(engine_.get());
    }
    int poisson(double mean) override { return CLHEP::RandPoisson::shoot(engine_.get(), mean); }

  private:
    void* enginePtr() override { return engine_.get(); }
    std::unique_ptr<CLHEP::HepRandomEngine> engine_;
  };
}  // namespace cepgen::clhep
using CLHEPRandomGenerator = cepgen::clhep::RandomGenerator;
REGISTER_RANDOM_GENERATOR("clhep", CLHEPRandomGenerator);
