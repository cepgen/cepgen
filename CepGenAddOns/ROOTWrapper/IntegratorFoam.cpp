#include "CepGen/Integration/Integrator.h"
#include "CepGen/Integration/Integrand.h"
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
  /// Integrand interfacing object used in Foam
  class FoamDensity : public TFoamIntegrand
  {
    public:
      explicit FoamDensity() = default;
      /// Specify the integrand
      void setIntegrand( Integrand* integr ) {
        integr_ = integr;
      }
      /// Compute the weight for a given phase space point
      inline double Density( int ndim, double* x ) override {
        if ( !integr_ )
          throw CG_FATAL( "FoamDensity" )
            << "Integrand object not yet initialised!";
        return integr_->eval( std::vector<double>( x, x+ndim ) );
      }

    private:
      Integrand* integr_; ///< Integrand
  };

  /// Foam general-purpose integration algorithm
  /// as developed by S. Jadach (Institute of Nuclear Physics, Krakow, PL)
  class IntegratorFoam : public Integrator
  {
    public:
      explicit IntegratorFoam( const ParametersList& );
      static std::string description() { return "FOAM general purpose MC integrator"; }

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
      foam_->SetkDim( integrand_->size() );
      density_->setIntegrand( integrand_ );
      foam_->ResetRho( density_.get() );
    CG_WARNING("")<<1;
      foam_->Initialize();
    CG_WARNING("")<<2;
      initialised_ = true;
    }
    for ( size_t i = 0; i < 100000; ++i )
      foam_->MakeEvent();
    //--- launch integration
    double norm, err;
    foam_->Finalize( norm, err );

    foam_->GetIntegMC( result, abserr );
    result_ = result;
    err_result_ = abserr;

    CG_DEBUG( "IntegratorFoam" ).log( [&]( auto& log ) {
      double eps = 5.e-4, avewt, wtmax, sigma;
      foam_->GetWtParams( eps, avewt, wtmax, sigma );
      const double ncalls = foam_->GetnCalls();
      const double effic = wtmax > 0 ? avewt/wtmax : 0.;
      log
        << "Result: " << result_ << " +- " << err_result_ << "\n\t"
        << "Relative error: " << err_result_/result_*100. << "%\n\t"
        << "Dispersion/<wt>= " << sigma
        << ", <wt>= " << avewt
        << ", <wt>/wtmax= " << effic << ",\n\t"
        << " for epsilon = " << eps << "\n\t"
        << " nCalls (initialisation only)= " << ncalls << ".";
    } );
  }
}

REGISTER_INTEGRATOR( "Foam", IntegratorFoam )
