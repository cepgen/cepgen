#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"
#include "cuba.h"

namespace cepgen {
  /// Cuba implementation of the Divonne integration algorithm
  class IntegratorCubaDivonne : public IntegratorCuba {
  public:
    explicit IntegratorCubaDivonne(const ParametersList &);
    void integrate(double &, double &) override;

    static ParametersDescription description();

  private:
    int key1_, key2_, key3_;
    int maxpass_;
    double border_, maxchisq_, mindeviation_;
    int ngiven_, ldxgiven_;
    int nextra_;
  };

  IntegratorCubaDivonne::IntegratorCubaDivonne(const ParametersList &params)
      : IntegratorCuba(params),
        key1_(steer<int>("KEY1")),
        key2_(steer<int>("KEY2")),
        key3_(steer<int>("KEY3")),
        maxpass_(steer<int>("MAXPASS")),
        border_(steer<double>("BORDER")),
        maxchisq_(steer<double>("MAXCHISQ")),
        mindeviation_(steer<double>("MINDEVIATION")),
        ngiven_(steer<int>("NGIVEN")),
        ldxgiven_(steer<int>("LDXGIVEN")),
        nextra_(steer<int>("NEXTRA")) {
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "Cuba-Divonne integrator built.";
  }

  void IntegratorCubaDivonne::integrate(double &result, double &abserr) {
    gIntegrand = integrand_;

    int nregions, neval, fail;
    double integral, error, prob;

    Divonne(integrand_->size(),
            1,
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
            ngiven_,
            ldxgiven_,
            nullptr,  // cubareal xgiven[]
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
    desc.add<int>("KEY1", 47);
    desc.add<int>("KEY2", 1);
    desc.add<int>("KEY3", 1);
    desc.add<int>("MAXPASS", 5);
    desc.add<double>("BORDER", 0.);
    desc.add<double>("MAXCHISQ", 10.);
    desc.add<double>("MINDEVIATION", 0.25);
    desc.add<int>("NGIVEN", 0);
    desc.add<int>("LDXGIVEN", 0);
    desc.add<int>("NEXTRA", 0);
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-divonne", IntegratorCubaDivonne)
