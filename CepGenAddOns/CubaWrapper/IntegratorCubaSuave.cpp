#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"
#include "cuba.h"

namespace cepgen {
  /// Cuba implementation of the VEGAS integration algorithm
  class IntegratorCubaSuave : public IntegratorCuba {
  public:
    explicit IntegratorCubaSuave(const ParametersList&);
    void integrate(double&, double&) override;

    static ParametersDescription description();

  private:
    int nnew_, nmin_;
    double flatness_;
    int verbose_;
  };

  IntegratorCubaSuave::IntegratorCubaSuave(const ParametersList& params)
      : IntegratorCuba(params),
        nnew_(steer<int>("NNEW")),
        nmin_(steer<int>("NMIN")),
        flatness_(steer<double>("FLATNESS")) {
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "Cuba-Suave integrator built.";
  }

  void IntegratorCubaSuave::integrate(double& result, double& abserr) {
    gIntegrand = integrand_;

    int neval, fail, nregions;
    double integral, error, prob;

    Suave(integrand_->size(),
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
          nnew_,
          nmin_,
          flatness_,
          nullptr,  // const char* statefile
          nullptr,  // void* spin
          &nregions,
          &neval,
          &fail,
          &integral,
          &error,
          &prob);

    result_ = result = integral;
    err_result_ = abserr = error;
  }

  ParametersDescription IntegratorCubaSuave::description() {
    auto desc = IntegratorCuba::description();
    desc.setDescription("Cuba implementation of the Suave algorithm");
    desc.add<int>("NNEW", 1000);
    desc.add<int>("NMIN", 2);
    desc.add<double>("FLATNESS", 25.);
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-suave", IntegratorCubaSuave)
