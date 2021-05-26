#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"

#include "cuba.h"

namespace cepgen {
  /// Cuba implementation of the VEGAS integration algorithm
  class IntegratorCubaSuave : public IntegratorCuba {
  public:
    explicit IntegratorCubaSuave(const ParametersList&);
    static std::string description() { return "Cuba implementation of the Suave algorithm"; }

    void integrate(double&, double&) override;

  private:
    int nnew_, nmin_;
    double flatness_;
    int verbose_;
  };

  IntegratorCubaSuave::IntegratorCubaSuave(const ParametersList& params)
      : IntegratorCuba(params),
        nnew_(params.get<int>("NNEW", 1000)),
        nmin_(params.get<int>("NMIN", 2)),
        flatness_(params.get<double>("FLATNESS", 25.)),
        verbose_(params.get<int>("verbose", 1)) {
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
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-suave", IntegratorCubaSuave)
