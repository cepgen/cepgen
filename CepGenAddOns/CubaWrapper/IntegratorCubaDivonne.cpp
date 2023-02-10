/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2022  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"
#include "cuba.h"

namespace cepgen {
  /// Cuba implementation of the Divonne integration algorithm
  class IntegratorCubaDivonne : public IntegratorCuba {
  public:
    explicit IntegratorCubaDivonne(const ParametersList&);
    void integrate(double&, double&) override;

    static ParametersDescription description();

  private:
    int key1_, key2_, key3_;
    int maxpass_;
    double border_, maxchisq_, mindeviation_;
    std::vector<std::vector<double> > given_;
    int ldxgiven_;
    int nextra_;
  };

  IntegratorCubaDivonne::IntegratorCubaDivonne(const ParametersList& params)
      : IntegratorCuba(params),
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
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "Cuba-Divonne integrator built.";
  }

  void IntegratorCubaDivonne::integrate(double& result, double& abserr) {
    int nregions, neval, fail;
    double integral, error, prob;
    int ngiven = given_.size();
    std::vector<double*> given_arr;
    for (auto& point : given_)
      given_arr.emplace_back(point.data());

    Divonne(integrand_->size(),
            ncomp_,
            cuba_integrand,
            nullptr,
            nvec_,
            epsrel_,
            epsabs_,
            verbose_,
            seed_,
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

    result_ = result = integral;
    err_result_ = abserr = error;
  }

  ParametersDescription IntegratorCubaDivonne::description() {
    auto desc = IntegratorCuba::description();
    desc.setDescription("Cuba implementation of the Divonne algorithm");
    desc.add<int>("Key1", 47).setDescription("sampling rule in the partitioning phase");
    desc.add<int>("Key2", 1).setDescription("sampling rule in the final integration phase");
    desc.add<int>("Key3", 1).setDescription(
        "strategy for the refinement phase"
        "(0 = do not treat the subregion any further, 1 = split the subregion up once more)");
    desc.add<int>("MaxPass", 5).setDescription("thoroughness parameter of the partitioning phase");
    desc.add<double>("Border", 0.).setDescription("border width of the integration region");
    desc.add<double>("MaxChisq", 10.)
        .setDescription(
            "maximum chi-square value a single subregion is allowed to have in the final integration phase");
    desc.add<double>("MinDeviation", 0.25)
        .setDescription(
            "fraction of the requested error of the entire integral, which determines whether it is worthwhile further "
            "examining a region that failed the chi-square test");
    desc.add<std::vector<std::vector<double> > >("Given", {})
        .setDescription("list of points where the integrand might have peaks");
    desc.add<int>("LDXGiven", 0)
        .setDescription("leading dimension of xgiven, i.e. the offset between one point and the next in memory");
    desc.add<int>("NExtra", 0).setDescription("maximum number of extra points the peak-finder subroutine will return");
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-divonne", IntegratorCubaDivonne);
