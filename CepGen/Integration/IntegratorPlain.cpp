#include "CepGen/Integration/IntegratorGSL.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"

#include <gsl/gsl_monte_plain.h>

namespace cepgen {
  /// Plain integration algorithm randomly sampling points in the phase space
  class IntegratorPlain : public IntegratorGSL {
  public:
    explicit IntegratorPlain(const ParametersList& params);
    static std::string description() { return "Plain (trial/error) integrator"; }

    void integrate(double&, double&) override;

  private:
    int ncvg_;
  };

  IntegratorPlain::IntegratorPlain(const ParametersList& params)
      : IntegratorGSL(params), ncvg_(params.get<int>("numFunctionCalls", 50000)) {}

  void IntegratorPlain::integrate(double& result, double& abserr) {
    //--- integration bounds
    std::vector<double> x_low(function_->dim, 0.), x_up(function_->dim, 1.);

    //--- launch integration
    std::unique_ptr<gsl_monte_plain_state, void (*)(gsl_monte_plain_state*)> pln_state(
        gsl_monte_plain_alloc(function_->dim), gsl_monte_plain_free);
    int res = gsl_monte_plain_integrate(
        function_.get(), &x_low[0], &x_up[0], function_->dim, ncvg_, rng_.get(), pln_state.get(), &result, &abserr);

    if (res != GSL_SUCCESS)
      throw CG_FATAL("Integrator:integrate") << "Error while performing the integration!\n\t"
                                             << "GSL error: " << gsl_strerror(res) << ".";

    result_ = result;
    err_result_ = abserr;
  }
}  // namespace cepgen

REGISTER_INTEGRATOR("plain", IntegratorPlain)
