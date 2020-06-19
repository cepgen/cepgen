#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"

#include <boost/math/quadrature/naive_monte_carlo.hpp>

namespace cepgen
{
  /// Boost's Naive integration algorithm
  class IntegratorNaive : public Integrator
  {
    public:
      explicit IntegratorNaive( const ParametersList& );

      void integrate( double&, double& ) override;

    private:
      std::function<double(const std::vector<double>&)> funct_;
      typedef boost::math::quadrature::naive_monte_carlo<double,decltype(funct_)> nmc_t;
      std::vector<std::pair<double,double> > bounds_;
      std::unique_ptr<nmc_t> mc_;
  };

  IntegratorNaive::IntegratorNaive( const ParametersList& params ) :
    Integrator( params ),
    funct_( [=]( const std::vector<double>& x ) -> double { return integrand_->eval( x ); } )
  {
    //--- a bit of printout for debugging
    CG_DEBUG( "Integrator:build" ) << "Boost's Naive integrator built.";
  }

  void
  IntegratorNaive::integrate( double& result, double& abserr )
  {
    if ( !initialised_ ) {
      bounds_ = std::vector<std::pair<double,double> >( function_->dim, { 0., 1. } );
      mc_.reset( new nmc_t( funct_, bounds_, 1.e-2, true, 1 ) );
      initialised_ = true;
    }

    auto task = mc_->integrate();

    result_ = result = task.get();
    err_result_ = abserr = mc_->current_error_estimate();
  }
}

REGISTER_INTEGRATOR( "Naive", IntegratorNaive )
