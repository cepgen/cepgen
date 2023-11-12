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
#include "CepGen/Utils/RandomGenerator.h"

namespace cepgen {
  class STLRandomGenerator : public utils::RandomGenerator {
  public:
    explicit STLRandomGenerator(const ParametersList& params) : utils::RandomGenerator(params) {
      const auto& type = steer<std::string>("type");
      if (type == "mt19937")
        rng_.reset(new Generator<std::mt19937>(rd_()));
      if (!rng_)
        throw CG_FATAL("STLRandomGenerator") << "Random number generator engine not set!";

      //gsl_rng_set(rng_.get(), seed_);

      CG_DEBUG("STLRandomGenerator") << "Random numbers generator: " << type << ".\n\t"
                                     << "Seed: " << seed_ << ".";
    }

    int uniformInt(int min, int max) override { return std::uniform_int_distribution<>(min, max)(*rng_); }

    double uniform(double min, double max) override { return std::uniform_real_distribution<>(min, max)(*rng_); }

    double normal(double mean, double rms) override { return std::normal_distribution<>(mean, rms)(*rng_); }

    double exponential(double exponent) override { return std::exponential_distribution<>(exponent)(*rng_); }

  private:
    std::random_device rd_;
    struct GeneratorObject {
      typedef std::random_device::result_type result_type;
      //long double operator()() { return 0.; }
      virtual long double min() = 0;
      virtual long double max() = 0;
    };
    template <typename T>
    class Generator : public GeneratorObject {
    public:
      explicit Generator(std::random_device::result_type rd) : gen_(rd) {}
      T& operator()() { return gen_; }
      long double min() override { return gen_.min(); }
      long double max() override { return gen_.max(); }

    private:
      T gen_;
    };
    std::unique_ptr<GeneratorObject> rng_;
  };
}  // namespace cepgen
