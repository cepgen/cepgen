#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"

#include "cuba.h"

namespace cepgen
{
  /// Boost naive integration algorithm
  class IntegratorCuba : public Integrator
  {
    public:
      explicit IntegratorCuba( const ParametersList& );

      void integrate( double&, double& ) override;

    private:
      int nvec_;
      double epsrel_, epsabs_;
      int mineval_, maxeval_;
      int nstart_, nincrease_, nbatch_;
      int gridno_;
      int verbose_;
  };

  IntegratorCuba::IntegratorCuba( const ParametersList& params ) :
    Integrator( params ),
    nvec_( params.get<int>( "NVEC", 1 ) ),
    epsrel_( params.get<double>( "EPSREL", 1.e-3 ) ),
    epsabs_( params.get<double>( "EPSABS", 1.e-12 ) ),
    mineval_( params.get<int>( "MINEVAL", 0 ) ),
    maxeval_( params.get<int>( "MAXEVAL", 50000 ) ),
    nstart_( params.get<int>( "NSTART", 1000 ) ),
    nincrease_( params.get<int>( "NINCREASE", 500 ) ),
    nbatch_( params.get<int>( "NBATCH", 1000 ) ),
    gridno_( params.get<int>( "GRIDNO", 0 ) ),
    verbose_( params.get<int>( "verbose", 1 ) )
  {
    //--- a bit of printout for debugging
    CG_DEBUG( "Integrator:build" ) << "CUBA integrator built.";
  }

  void
  IntegratorCuba::integrate( double& result, double& abserr )
  {
    auto integr = [&]( const int* ndim, const cubareal xx[], const int* ncomp, cubareal ff[], void* userdata ) -> int {
      ff[0] = function_->f( (double*)xx, function_->dim, (void*)function_->params );
      return 0;
    };

    int neval, fail;
    cubareal integral, error, prob;

    Vegas( function_->dim, 1, integr, nullptr, nvec_,
      epsrel_, epsabs_, verbose_, seed_,
      mineval_, maxeval_, nstart_, nincrease_, nbatch_,
      gridno_, nullptr, nullptr,
      &neval, &fail, &integral, &error, &prob);

    result_ = result = integral;
    err_result_ = abserr = error;
  }
}

REGISTER_INTEGRATOR( "cuba", IntegratorCuba )
