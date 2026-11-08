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

#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
#include <TRandomGen.h>

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace std::string_literals;

namespace cepgen::root {
  class RandomGenerator final : public utils::RandomGenerator {
  public:
    explicit RandomGenerator(const ParametersList& params) : utils::RandomGenerator(params) {
      if (const auto& type = steer<std::string>("type"); type == "Ranlux")
        random_number_generator_ = std::make_unique<TRandom1>();
      else if (type == "Tausworthe")
        random_number_generator_ = std::make_unique<TRandom2>();
      else if (type == "MersenneTwister")
        random_number_generator_ = std::make_unique<TRandom3>();
      else if (type == "Ranluxpp")
        random_number_generator_ = std::make_unique<TRandomRanluxpp>();
      else if (type == "MixMax")
        random_number_generator_ = std::make_unique<TRandomMixMax>();
      else if (type == "MixMax17")
        random_number_generator_ = std::make_unique<TRandomMixMax17>();
      else if (type == "MixMax256")
        random_number_generator_ = std::make_unique<TRandomMixMax256>();
      else
        throw CG_FATAL("root:RandomGenerator") << "Random number generator engine invalid: '" << type << "'.";
      random_number_generator_->SetSeed(seed_);
    }

    static ParametersDescription description() {
      auto desc = utils::RandomGenerator::description();
      desc.setDescription("ROOT random number generator engine");
      desc.add("type", "Ranlux"s)
          .allow("Ranlux")
          .allow("Tausworthe")
          .allow("MersenneTwister")
          .allow("Ranluxpp")
          .allow("MixMax")
          .allow("MixMax17")
          .allow("MixMax256")
          .setDescription("random number engine");
      return desc;
    }

    int uniformInt(int min, int max) override { return min + random_number_generator_->Integer(max - min + 1); }
    double uniform(double min, double max) override { return random_number_generator_->Uniform(min, max); }
    double normal(double mean, double rms) override { return random_number_generator_->Gaus(mean, rms); }
    double exponential(double exponent) override { return random_number_generator_->Exp(exponent); }
    double breitWigner(double mean, double scale) override {
      return random_number_generator_->BreitWigner(mean, scale);
    }
    double landau(double location, double width) override { return random_number_generator_->Landau(location, width); }
    int poisson(double mean) override { return random_number_generator_->Poisson(mean); }

  private:
    void* enginePtr() override { return random_number_generator_.get(); }
    std::unique_ptr<TRandom> random_number_generator_;
  };
}  // namespace cepgen::root
using ROOTRandomGenerator = cepgen::root::RandomGenerator;
REGISTER_RANDOM_GENERATOR("root", ROOTRandomGenerator);
