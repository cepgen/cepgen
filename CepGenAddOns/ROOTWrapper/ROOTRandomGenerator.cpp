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

#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandomGen.h>

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen {
  class ROOTRandomGenerator : public utils::RandomGenerator {
  public:
    explicit ROOTRandomGenerator(const ParametersList& params) : utils::RandomGenerator(params) {
      const auto& type = steer<std::string>("type");
      if (type == "Ranlux")
        rng_.reset(new TRandom1);
      else if (type == "Tausworthe")
        rng_.reset(new TRandom2);
      else if (type == "MersenneTwister")
        rng_.reset(new TRandom3);
      else if (type == "Ranlux")
        rng_.reset(new TRandomRanluxpp);
      else if (type == "MixMax")
        rng_.reset(new TRandomMixMax);
      else if (type == "MixMax17")
        rng_.reset(new TRandomMixMax17);
      else if (type == "MixMax256")
        rng_.reset(new TRandomMixMax256);
      else
        throw CG_FATAL("ROOTRandomGenerator") << "Random number generator engine invalid: '" << type << "'.";

      rng_->SetSeed(seed_);
    }

    static ParametersDescription description() {
      auto desc = utils::RandomGenerator::description();
      desc.setDescription("ROOT random number generator engine");
      desc.add<std::string>("type", "Ranlux").setDescription("random number engine");
      return desc;
    }

    int uniformInt(int min, int max) override { return min + rng_->Integer(max - min + 1); }
    double uniform(double min, double max) override { return rng_->Uniform(min, max); }
    double normal(double mean, double rms) override { return rng_->Gaus(mean, rms); }
    double exponential(double exponent) override { return rng_->Exp(exponent); }
    double breitWigner(double mean, double scale) override { return rng_->BreitWigner(mean, scale); }
    double landau(double location, double width) override { return rng_->Landau(location, width); }
    int poisson(double mean) override { return rng_->Poisson(mean); }

  private:
    void* enginePtr() override { return rng_.get(); }
    std::unique_ptr<TRandom> rng_;
  };
}  // namespace cepgen

REGISTER_RANDOM_GENERATOR("root", ROOTRandomGenerator);
