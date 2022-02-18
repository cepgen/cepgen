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
        nnew_(steer<int>("NNew")),
        nmin_(steer<int>("NMin")),
        flatness_(steer<double>("Flatness")) {
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
    desc.add<int>("NNew", 1000).setDescription("number of new integrand evaluations in each subdivision");
    desc.add<int>("NMin", 2).setDescription(
        "minimum number of samples a former pass must contribute to a subregion to be considered in that region’s "
        "compound integral value");
    desc.add<double>("Flatness", 50.).setDescription("type of norm used to compute the fluctuation of a sample");
    return desc;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-suave", IntegratorCubaSuave)
