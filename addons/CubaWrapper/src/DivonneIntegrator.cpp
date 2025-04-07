/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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

#include <cuba.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGenCuba/Integrator.h"

namespace cepgen::cuba {
  /// Cuba implementation of the Divonne integration algorithm
  class DivonneIntegrator : public Integrator {
  public:
    explicit DivonneIntegrator(const ParametersList& params)
        : Integrator(params),
          key1_(steer<int>("Key1")),
          key2_(steer<int>("Key2")),
          key3_(steer<int>("Key3")),
          maxpass_(steer<int>("MaxPass")),
          border_(steer<double>("Border")),
          maxchisq_(steer<double>("MaxChisq")),
          mindeviation_(steer<double>("MinDeviation")),
          given_(steer<std::vector<std::vector<double> > >("Given")),
          ldxgiven_(steer<int>("LDXGiven")),
          nextra_(steer<int>("NExtra")) {
      CG_DEBUG("Integrator:build") << "Cuba-Divonne integrator built.";
    }

    static ParametersDescription description() {
      auto desc = Integrator::description();
      desc.setDescription("Cuba implementation of the Divonne algorithm");
      desc.add("Key1", 47).setDescription("sampling rule in the partitioning phase");
      desc.add("Key2", 1).setDescription("sampling rule in the final integration phase");
      desc.add("Key3", 1)
          .allow(0, "do not treat the subregion any further")
          .allow(1, "split the subregion up once more")
          .setDescription("strategy for the refinement phase");
      desc.add("MaxPass", 5).setDescription("thoroughness parameter of the partitioning phase");
      desc.add("Border", 0.).setDescription("border width of the integration region");
      desc.add("MaxChisq", 10.)
          .setDescription(
              "maximum chi-square value a single subregion is allowed to have in the final integration phase");
      desc.add("MinDeviation", 0.25)
          .setDescription(
              "fraction of the requested error of the entire integral, which determines whether it is worthwhile "
              "further examining a region that failed the chi-square test");
      desc.add("Given", std::vector<std::vector<double> >{})
          .setDescription("list of points where the integrand might have peaks");
      desc.add("LDXGiven", 0)
          .setDescription("leading dimension of xgiven, i.e. offset between one point and next in memory");
      desc.add("NExtra", 0).setDescription("maximum number of extra points the peak-finder subroutine will return");
      return desc;
    }

    Value integrate() override {
      int nregions, neval, fail;
      double integral, error, prob;
      int ngiven = given_.size();
      std::vector<double*> given_arr;
      std::transform(
          given_.begin(), given_.end(), std::back_inserter(given_arr), [](auto& point) { return point.data(); });

      Divonne(gIntegrand->size(),
              ncomp_,
              cuba_integrand,
              nullptr,
              nvec_,
              epsrel_,
              epsabs_,
              verbosity_,
              steerAs<unsigned long long, int>("seed"),
              mineval_,
              maxeval_,
              key1_,
              key2_,
              key3_,
              maxpass_,
              border_,
              maxchisq_,
              mindeviation_,
              ngiven,
              ldxgiven_,
              !given_arr.empty() ? *given_arr.data() : nullptr,  // cubareal xgiven[]
              nextra_,
              nullptr,  // peakfinder_t peakfinder
              nullptr,  // const char *statefile
              nullptr,  // void *spin
              &nregions,
              &neval,
              &fail,
              &integral,
              &error,
              &prob);

      return Value{integral, error};
    }

  private:
    const int key1_, key2_, key3_;
    const int maxpass_;
    const double border_, maxchisq_, mindeviation_;
    std::vector<std::vector<double> > given_;
    const int ldxgiven_;
    const int nextra_;
  };
}  // namespace cepgen::cuba
using cepgen::cuba::DivonneIntegrator;
REGISTER_INTEGRATOR("cuba_divonne", DivonneIntegrator);
