#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"

#include "cuba.h"

namespace cepgen {
  /// Cuba implementation of the Divonne integration algorithm
  class IntegratorCubaDivonne : public IntegratorCuba {
  public:
    explicit IntegratorCubaDivonne(const ParametersList &);
    static std::string description() { return "Cuba implementation of the Divonne algorithm"; }

    void integrate(double &, double &) override;

  private:
    int key1_, key2_, key3_;
    int maxpass_;
    double border_, maxchisq_, mindeviation_;
    int ngiven_, ldxgiven_;
    int nextra_;
  };

  IntegratorCubaDivonne::IntegratorCubaDivonne(const ParametersList &params)
      : IntegratorCuba(params),
        key1_(params.get<int>("KEY1", 47)),
        key2_(params.get<int>("KEY2", 1)),
        key3_(params.get<int>("KEY3", 1)),
        maxpass_(params.get<int>("MAXPASS", 5)),
        border_(params.get<double>("BORDER", 0.)),
        maxchisq_(params.get<double>("MAXCHISQ", 10.)),
        mindeviation_(params.get<double>("MINDEVIATION", 0.25)),
        ngiven_(params.get<int>("NGIVEN", 0)),
        ldxgiven_(params.get<int>("LDXGIVEN", 0)),
        nextra_(params.get<int>("NEXTRA", 0)) {
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
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-divonne", IntegratorCubaDivonne)
