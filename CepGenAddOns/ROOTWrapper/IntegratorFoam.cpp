#include "CepGen/Integration/Integrator.h"
#include "CepGen/Modules/IntegratorFactory.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"

#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"

#include "TFoam.h"
#include "TFoamIntegrand.h"

namespace cepgen
{
  /// Integrand interfacing object used in FOAM
  class FoamDensity : public TFoamIntegrand
  {
    public:
      explicit FoamDensity() = default;
      /// Specify the runtime parameters for the integrand steering
      void setParameters( Parameters* params ) {
        params_.reset( params );
      }
      /// Specify the integrand functional
      inline void setFunction( double(*func)( double*, unsigned long, void* ) ) {
        function_ = func;
      }
      /// Compute the weight for a given phase space point
      inline double Density( int ndim, double* x ) override {
        return function_( x, ndim, (void*)params_.get() );
      }

    private:
      std::shared_ptr<Parameters> params_; ///< Runtime parameters
      std::function<double( double*, size_t, void* )> function_; ///< Integrand
  };

  /// FOAM general-purpose integration algorithm
  class IntegratorFoam : public Integrator
  {
    public:
      explicit IntegratorFoam( const ParametersList& );

      void integrate( double&, double& ) override;
      inline double uniform() const override { return rnd_->Rndm(); }

    private:
      std::unique_ptr<TFoam> foam_;
      std::unique_ptr<TRandom> rnd_;
      std::unique_ptr<FoamDensity> density_;
  };

  IntegratorFoam::IntegratorFoam( const ParametersList& params ) :
    Integrator( params ),
    foam_( new TFoam( "Foam" ) ), density_( new FoamDensity )
  {
    const auto& rnd_mode = params.get<std::string>( "rngEngine", "MersenneTwister" );
    if ( rnd_mode == "Ranlux" )
      rnd_.reset( new TRandom1 );
    else if ( rnd_mode == "generic" )
      rnd_.reset( new TRandom2 );
    else if ( rnd_mode == "MersenneTwister" )
      rnd_.reset( new TRandom3 );
    else
      throw CG_FATAL( "IntegratorFoam" )
        << "Unrecognised random generator: \"" << rnd_mode << "\".";
    rnd_->SetSeed( seed_ );

    foam_->SetPseRan( rnd_.get() );
    foam_->SetnCells( params.get<int>( "nCells", 1000 ) );
    foam_->SetnSampl( params.get<int>( "nSampl", 200 ) );
    foam_->SetnBin( params.get<int>( "nBin", 8 ) );
    foam_->SetEvPerBin( params.get<int>( "EvPerBin", 25 ) );
    foam_->SetChat( verbosity_ );
    //--- a bit of printout for debugging
    CG_DEBUG( "Integrator:build" ) << "FOAM integrator built\n\t"
      << "Version: " << foam_->GetVersion() << ".";
  }

  void
  IntegratorFoam::integrate( double& result, double& abserr )
  {
    if ( !initialised_ ) {
      density_->setParameters( (Parameters*)function_->params );
      density_->setFunction( function_->f );
      foam_->SetkDim( function_->dim );
      foam_->SetRho( density_.get() );
      foam_->Initialize();
      initialised_ = true;
    }
    for ( size_t i = 0; i < 100000; ++i )
      foam_->MakeEvent();
    //--- launch integration
    foam_->GetIntegMC( result, abserr );
    double norm, err;
    foam_->Finalize( norm, err );

    result_ = result;
    err_result_ = abserr;
  }
}

REGISTER_INTEGRATOR( "Foam", IntegratorFoam )
