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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace std::string_literals;

class GSLRandomGenerator final : public utils::RandomGenerator {
public:
  explicit GSLRandomGenerator(const ParametersList& params) : RandomGenerator(params) {
    gsl_rng_env_setup();
    gsl_rng_type* rng_engine;
    const auto& type = steer<std::string>("type");
    if (type == "mt19937"s)
      rng_engine = const_cast<gsl_rng_type*>(gsl_rng_mt19937);
    else if (type == "taus"s)
      rng_engine = const_cast<gsl_rng_type*>(gsl_rng_taus);
    else if (type == "taus2"s)
      rng_engine = const_cast<gsl_rng_type*>(gsl_rng_taus2);
    else if (type == "gfsr4"s)
      rng_engine = const_cast<gsl_rng_type*>(gsl_rng_gfsr4);
    else if (type == "ranlxs0"s)
      rng_engine = const_cast<gsl_rng_type*>(gsl_rng_ranlxs0);
    else
      throw CG_FATAL("GSLRandomGenerator") << "Random number generator engine invalid: '" << type << "'.";

    rng_.reset(gsl_rng_alloc(rng_engine));
    gsl_rng_set(rng_.get(), seed_);

    CG_DEBUG("GSLRandomGenerator") << "Random numbers generator: " << gsl_rng_name(rng_.get()) << ". Seed: " << seed_
                                   << ".";
  }

  static ParametersDescription description() {
    auto desc = utils::RandomGenerator::description();
    desc.setDescription("GSL random number generator engine");
    desc.add("type", "mt19937"s)
        .allow("mt19937"s, "Mersenne-Twister generator")
        .allow("taus"s, "maximally equi-distributed combined Tausworthe generator by L’Écuyer")
        .allow("taus2"s,
               "maximally equi-distributed combined Tausworthe generator by L’Écuyer (w/ improved seeding procedure)")
        .allow("gfsr4"s, "lagged-fibonacci generator")
        .allow("ranlxs0"s, "second-generation version of the RANLUX algorithm of Luscher")
        .setDescription("random number engine");
    return desc;
  }

  int uniformInt(int min, int max) override { return min + gsl_rng_uniform_int(rng_.get(), max - min + 1); }
  double uniform(double min, double max) override { return Limits{min, max}.x(gsl_rng_uniform(rng_.get())); }
  double normal(double mean, double rms) override { return gsl_ran_gaussian(rng_.get(), rms) + mean; }
  double exponential(double exponent) override { return gsl_ran_exponential(rng_.get(), exponent); }
  double breitWigner(double mean, double scale) override { return gsl_ran_cauchy(rng_.get(), scale) + mean; }
  double landau(double location, double width) override { return width * gsl_ran_landau(rng_.get()) + location; }
  int poisson(double mean) override { return gsl_ran_poisson(rng_.get(), mean); }

private:
  /// A deleter object for GSL's random number generator
  struct gsl_rng_deleter {
    /// Destructor method for the random number generator service
    inline void operator()(gsl_rng* rng) const { gsl_rng_free(rng); }
  };
  std::unique_ptr<gsl_rng, gsl_rng_deleter> rng_;  ///< Instance of random number generator service

  void* enginePtr() override { return rng_.get(); }
};
REGISTER_RANDOM_GENERATOR("gsl", GSLRandomGenerator);
