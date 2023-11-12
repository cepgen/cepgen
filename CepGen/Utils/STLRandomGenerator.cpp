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
  template <typename T>
  class STLRandomGenerator : public utils::RandomGenerator {
  public:
    explicit STLRandomGenerator(const ParametersList& params) : utils::RandomGenerator(params) {
      std::random_device rd;
      rng_.reset(new T(seed_ > 0ull ? seed_ : rd()));
      if (!rng_)
        throw CG_FATAL("STLRandomGenerator") << "Random number generator engine not set!";

      CG_DEBUG("STLRandomGenerator") << "Random numbers generator with seed: " << seed_ << ".";
    }

    int uniformInt(int min, int max) override { return std::uniform_int_distribution<>(min, max)(*rng_); }

    double uniform(double min, double max) override { return std::uniform_real_distribution<>(min, max)(*rng_); }

    double normal(double mean, double rms) override { return std::normal_distribution<>(mean, rms)(*rng_); }

    double exponential(double exponent) override { return std::exponential_distribution<>(exponent)(*rng_); }

  private:
    std::unique_ptr<T> rng_;
  };
}  // namespace cepgen

typedef cepgen::STLRandomGenerator<std::mt19937> mt19937;
typedef cepgen::STLRandomGenerator<std::mt19937_64> mt19937_64;
typedef cepgen::STLRandomGenerator<std::ranlux24_base> ranlux24_base;
typedef cepgen::STLRandomGenerator<std::ranlux48_base> ranlux48_base;
REGISTER_RANDOM_GENERATOR("stl:mt19937", mt19937);
REGISTER_RANDOM_GENERATOR("stl:mt19937_64", mt19937_64);
REGISTER_RANDOM_GENERATOR("stl:ranlux24_base", ranlux24_base);
REGISTER_RANDOM_GENERATOR("stl:ranlux48_base", ranlux48_base);
