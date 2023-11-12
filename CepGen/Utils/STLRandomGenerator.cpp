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

#include <memory>
#include <random>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen {
  class STLRandomGenerator : public utils::RandomGenerator {
  public:
    explicit STLRandomGenerator(const ParametersList& params) : utils::RandomGenerator(params) {
      std::random_device rd;
      const auto seed = seed_ > 0ull ? seed_ : rd();
      const auto& type = steer<std::string>("type");
      if (type == "mt19937")
        gen_.reset(new Generator<std::mt19937>(seed));
      else if (type == "mt19937_64")
        gen_.reset(new Generator<std::mt19937_64>(seed));
      else if (type == "ranlux24_base")
        gen_.reset(new Generator<std::ranlux24_base>(seed));
      else if (type == "ranlux48_base")
        gen_.reset(new Generator<std::ranlux48_base>(seed));
      else
        throw CG_FATAL("STLRandomGenerator") << "Random number generator engine not set!";

      CG_DEBUG("STLRandomGenerator") << "Random numbers generator with seed: " << seed_ << ".";
    }

    static ParametersDescription description() {
      auto desc = utils::RandomGenerator::description();
      desc.setDescription("STL random number generator engine");
      desc.add<std::string>("type", "mt19937").setDescription("random number engine");
      return desc;
    }

    int uniformInt(int min, int max) override { return gen_->uniformInt(min, max); }
    double uniform(double min, double max) override { return gen_->uniform(min, max); }
    double normal(double mean, double rms) override { return gen_->normal(mean, rms); }
    double exponential(double exponent) override { return gen_->exponential(exponent); }

  private:
    template <typename T>
    class Generator : public utils::RandomGenerator {
    public:
      explicit Generator(unsigned long int value) : utils::RandomGenerator(ParametersList()), rng_(value) {}
      int uniformInt(int min, int max) override { return std::uniform_int_distribution<>(min, max)(rng_); }
      double uniform(double min, double max) override { return std::uniform_real_distribution<>(min, max)(rng_); }
      double normal(double mean, double rms) override { return std::normal_distribution<>(mean, rms)(rng_); }
      double exponential(double exponent) override { return std::exponential_distribution<>(exponent)(rng_); }

    private:
      T rng_;
    };

    std::unique_ptr<utils::RandomGenerator> gen_;
  };
}  // namespace cepgen

REGISTER_RANDOM_GENERATOR("stl", STLRandomGenerator);
