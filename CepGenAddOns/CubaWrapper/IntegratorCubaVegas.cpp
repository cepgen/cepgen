#include "CepGenAddOns/CubaWrapper/IntegratorCuba.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"

#include "cuba.h"

namespace cepgen {
  /// Cuba implementation of the VEGAS integration algorithm
  class IntegratorCubaVegas : public IntegratorCuba {
  public:
    explicit IntegratorCubaVegas(const ParametersList&);
    static std::string description() { return "Cuba implementation of the VEGAS algorithm"; }

    void integrate(double&, double&) override;

  private:
    int nstart_, nincrease_, nbatch_;
    int gridno_;
    int verbose_;
  };

  IntegratorCubaVegas::IntegratorCubaVegas(const ParametersList& params)
      : IntegratorCuba(params),
        nstart_(params.get<int>("NSTART", 1000)),
        nincrease_(params.get<int>("NINCREASE", 500)),
        nbatch_(params.get<int>("NBATCH", 1000)),
        gridno_(params.get<int>("GRIDNO", 0)),
        verbose_(params.get<int>("verbose", 1)) {
    //--- a bit of printout for debugging
    CG_DEBUG("Integrator:build") << "Cuba-VEGAS integrator built.";
  }

  void IntegratorCubaVegas::integrate(double& result, double& abserr) {
    gIntegrand = integrand_;
    /*auto integr = [&]( const int* ndim, const double xx[], const int* ncomp, double ff[], void* userdata ) -> int {
      //ff[0] = function_->f( (double*)xx, function_->dim, (void*)function_->params );
      ff[0] = integrand_->eval( std::vector<double>( xx, xx+*ndim ) );
      return 0;
    };*/

    int neval, fail;
    double integral, error, prob;

    Vegas(integrand_->size(),
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
          nstart_,
          nincrease_,
          nbatch_,
          gridno_,
          nullptr,
          nullptr,
          &neval,
          &fail,
          &integral,
          &error,
          &prob);

    result_ = result = integral;
    err_result_ = abserr = error;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("cuba-vegas", IntegratorCubaVegas)
