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

#include <memory>
#include <random>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace std::string_literals;

class STLRandomGenerator final : public utils::RandomGenerator {
public:
  explicit STLRandomGenerator(const ParametersList& params) : RandomGenerator(params) {
    std::random_device rd;
    const auto seed = seed_ > 0ull ? seed_ : rd();
    if (const auto& type = steer<std::string>("type"); type == "default"s)
      gen_ = std::make_unique<Generator<std::default_random_engine> >(seed);
    else if (type == "minstd_rand0"s)
      gen_ = std::make_unique<Generator<std::minstd_rand0> >(seed);
    else if (type == "minstd_rand"s)
      gen_ = std::make_unique<Generator<std::minstd_rand> >(seed);
    else if (type == "mt19937"s)
      gen_ = std::make_unique<Generator<std::mt19937> >(seed);
    else if (type == "mt19937_64"s)
      gen_ = std::make_unique<Generator<std::mt19937_64> >(seed);
    else if (type == "ranlux24_base"s)
      gen_ = std::make_unique<Generator<std::ranlux24_base> >(seed);
    else if (type == "ranlux48_base"s)
      gen_ = std::make_unique<Generator<std::ranlux48_base> >(seed);
    else if (type == "ranlux24"s)
      gen_ = std::make_unique<Generator<std::ranlux24> >(seed);
    else if (type == "ranlux48"s)
      gen_ = std::make_unique<Generator<std::ranlux48> >(seed);
    else if (type == "knuth_b"s)
      gen_ = std::make_unique<Generator<std::knuth_b> >(seed);
    else
      throw CG_FATAL("STLRandomGenerator") << "Random number generator engine invalid: '" << type << "'.";

    CG_DEBUG("STLRandomGenerator") << "Random numbers generator with seed: " << seed_ << ".";
  }

  static ParametersDescription description() {
    auto desc = RandomGenerator::description();
    desc.setDescription("STL random number generator engine");
    desc.add("type", "default"s)
        .allow("default", "implementation-defined algorithm")
        .allow("minstd_rand0",
               "Discovered in 1969 by Lewis, Goodman and Miller, adopted as \"Minimal standard\" in 1988 by Park and "
               "Miller")
        .allow("minstd_rand", "Newer \"Minimum standard\", recommended by Park, Miller, and Stockmeyer in 1993")
        .allow("mt19937", "32-bit Mersenne Twister by Matsumoto and Nishimura, 1998")
        .allow("mt19937_64", "64-bit Mersenne Twister by Matsumoto and Nishimura, 2000")
        .allow("ranlux24_base", "subtract-w/-carry algorithm (24, 10, 24)")
        .allow("ranlux48_base", "subtract-w/-carry algorithm (48, 5, 12)")
        .allow("ranlux24", "24-bit RANLUX generator by Martin Lüscher and Fred James, 1994")
        .allow("ranlux48", "48-bit RANLUX generator by Martin Lüscher and Fred James, 1994")
        .allow("knuth_b", "PRN engine adaptor discarding a certain amount of data produced by base engine (389, 11)")
        .setDescription("random number engine");
    return desc;
  }

  int uniformInt(int min, int max) override { return gen_->uniformInt(min, max); }
  double uniform(double min, double max) override { return gen_->uniform(min, max); }
  double normal(double mean, double rms) override { return gen_->normal(mean, rms); }
  double exponential(double exponent) override { return gen_->exponential(exponent); }
  double breitWigner(double mean, double scale) override { return gen_->breitWigner(mean, scale); }
  int poisson(double mean) override { return gen_->poisson(mean); }

private:
  template <typename T>
  class Generator : public RandomGenerator {
  public:
    explicit Generator(unsigned long int value) : RandomGenerator(ParametersList()), rng_(value) {}
    int uniformInt(int min, int max) override { return std::uniform_int_distribution<>(min, max)(rng_); }
    double uniform(double min, double max) override { return std::uniform_real_distribution<>(min, max)(rng_); }
    double normal(double mean, double rms) override { return std::normal_distribution<>(mean, rms)(rng_); }
    double exponential(double exponent) override { return std::exponential_distribution<>(exponent)(rng_); }
    double breitWigner(double mean, double scale) override { return std::cauchy_distribution<>(mean, scale)(rng_); }
    int poisson(double mean) override { return std::poisson_distribution<>(mean)(rng_); }

  private:
    T rng_;
  };

  std::unique_ptr<RandomGenerator> gen_;
};
REGISTER_RANDOM_GENERATOR("stl", STLRandomGenerator);
