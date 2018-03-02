#include "GenericKTProcess.h"
#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace Process
  {
    GenericKTProcess::GenericKTProcess( const std::string& name,
                                        const std::string& description,
                                        const unsigned int& num_user_dimensions,
                                        const std::array<ParticleCode,2>& partons,
                                        const std::vector<ParticleCode>& central ) :
      GenericProcess( name, description+" (kT-factorisation approach)" ),
      kt_jacobian_( 0. ), central_jacobian_( 0. ),
      log_qmin_( 0. ), log_qmax_( 0. ),
      qt1_( 0. ), phi_qt1_( 0. ),
      qt2_( 0. ), phi_qt2_( 0. ),
      flux1_( 0. ), flux2_( 0. ),
      kNumUserDimensions( num_user_dimensions ),
      kIntermediateParts( partons ), kProducedParts( central )
    {}

    void
    GenericKTProcess::prepareKTKinematics()
    {
      DebuggingInsideLoop("Dummy kinematics prepared!");
    }

    void
    GenericKTProcess::addEventContent()
    {
      GenericProcess::setEventContent(
        { // incoming state
          { Particle::IncomingBeam1, Proton },
          { Particle::IncomingBeam2, Proton },
          { Particle::Parton1, kIntermediateParts[0] },
          { Particle::Parton2, kIntermediateParts[1] }
        },
        { // outgoing state
          { Particle::OutgoingBeam1, { Proton } },
          { Particle::OutgoingBeam2, { Proton } },
          { Particle::CentralSystem, kProducedParts }
        }
      );
      setExtraContent();
    }

    unsigned int
    GenericKTProcess::numDimensions( const Kinematics::ProcessMode& process_mode ) const
    {
      switch ( process_mode ) {
        default:
        case Kinematics::ElasticElastic:     return kNumRequiredDimensions+kNumUserDimensions;
        case Kinematics::ElasticInelastic:
        case Kinematics::InelasticElastic:   return kNumRequiredDimensions+kNumUserDimensions+1;
        case Kinematics::InelasticInelastic: return kNumRequiredDimensions+kNumUserDimensions+2;
      }
    }

    void
    GenericKTProcess::addPartonContent()
    {
      // Incoming partons
      qt1_ = exp( log_qt_limits_.x( x( 0 ) ) );
      qt2_ = exp( log_qt_limits_.x( x( 1 ) ) );
      phi_qt1_ = phi_qt_limits_.x( x( 2 ) );
      phi_qt2_ = phi_qt_limits_.x( x( 3 ) );
      DebuggingInsideLoop( Form( "photons transverse virtualities (qt):\n\t"
                                 "  mag = %g / %g (%g < log(qt) < %g)\n\t"
                                 "  phi = %g / %g (%g < Dphi < %g)",
                                 qt1_, qt2_, log_qt_limits_.min(), log_qt_limits_.max(),
                                 phi_qt1_, phi_qt2_, phi_qt_limits_.min(), phi_qt_limits_.max() ) );
    }

    void
    GenericKTProcess::setKinematics( const Kinematics& kin )
    {
      cuts_ = kin;

      //----------------------------------------------------------------

      kt_jacobian_ = 1.;

      //----------------------------------------------------------------

      if ( cuts_.cuts.initial.count( Cuts::qt ) == 0
        || !cuts_.cuts.initial.at( Cuts::qt ).valid() ) {
        InWarning( "Failed to retrieve a transverse virtuality range for the incoming partons from the user configuration!\n\t"
                   "Setting it to the default 10¯¹⁰ < qT < 500 GeV value." );
        cuts_.cuts.initial[Cuts::qt] = { 1.e-10, 500. };
      }
      Kinematics::Limits qt_limits = cuts_.cuts.initial.at( Cuts::qt );
      log_qt_limits_ = {
        std::max( log( qt_limits.min() ), -10. ), //FIXME fine-tuning needed!
        log( qt_limits.max() )
      };

      if ( log_qt_limits_.max() > 10. )
        InWarning( Form( "qT range above \"reasonable\" range: maximal qT specified = 10^(%.1f)", log_qt_limits_.max() ) );

      kt_jacobian_ *= pow( log_qt_limits_.range(), 2 ); // d(q1t) d(q2t)

      //----------------------------------------------------------------

      if ( cuts_.cuts.initial.count( Cuts::phi_qt ) == 0
        || !cuts_.cuts.initial.at( Cuts::phi_qt ).valid() ) {
        InWarning( "Failed to retrieve an incoming partons azimuthal angle difference from the user configuration!\n\t"
                   "Setting it to the default 0 < Δɸ < 2π value." );
        cuts_.cuts.initial[Cuts::phi_qt] = { 0., 2.*M_PI };
      }
      phi_qt_limits_ = cuts_.cuts.initial.at( Cuts::phi_qt );
      kt_jacobian_ *= pow( phi_qt_limits_.range(), 2 ); // d(phi1) d(phi2)

      //----------------------------------------------------------------

      if ( cuts_.cuts.remnants.count( Cuts::mass ) == 0
        || !cuts_.cuts.remnants.at( Cuts::mass ).valid() ) {
        InWarning( "Failed to retrieve a diffractive system mass range from the user configuration!\n\t"
                   "Setting it to the default 1.07 < Mx < 1000 GeV value." );
        cuts_.cuts.remnants[Cuts::mass] = { 1.07, 1000. };
      }
      mx_limits_ = cuts_.cuts.remnants.at( Cuts::mass );
      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic: default:
          break;
        case Kinematics::ElasticInelastic:
          kt_jacobian_ *= 2.* mx_limits_.range();
          break;
        case Kinematics::InelasticElastic:
          kt_jacobian_ *= 2.* mx_limits_.range();
          break;
        case Kinematics::InelasticInelastic:
          kt_jacobian_ *= 2.* mx_limits_.range();
          kt_jacobian_ *= 2.* mx_limits_.range();
          break;
      } // d(mx/y**2)

      preparePhaseSpace();
    }

    double
    GenericKTProcess::computeWeight()
    {
      addPartonContent();
      computeOutgoingPrimaryParticlesMasses();
      prepareKTKinematics();

      double res_jacobian = qt1_*qt2_;
      switch ( cuts_.mode ) {
        case Kinematics::InelasticElastic:   res_jacobian = MX_; break;
        case Kinematics::ElasticInelastic:   res_jacobian = MY_; break;
        case Kinematics::InelasticInelastic: res_jacobian = MX_*MY_; break;
        default: break;
      }

      const double integrand = computeKTFactorisedMatrixElement();
      const double weight = ( kt_jacobian_*central_jacobian_*res_jacobian ) * integrand;

      DebuggingInsideLoop( Form( "\n\tJacobian = %g * %g * %g\n\tIntegrand = %g\n\tdW = %g", kt_jacobian_, central_jacobian_, res_jacobian, integrand, weight ) );

      return weight;
    }

    void
    GenericKTProcess::computeOutgoingPrimaryParticlesMasses()
    {
      const unsigned int op_index = kNumRequiredDimensions+kNumUserDimensions;
      switch ( cuts_.mode ) {
        case Kinematics::ElectronProton: default: {
          InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
          exit( 0 ); } break;
        case Kinematics::ElasticElastic: {
          MX_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
          MY_ = event_->getOneByRole( Particle::IncomingBeam2 ).mass();
        } break;
        case Kinematics::ElasticInelastic: {
          MX_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
          MY_ = mx_limits_.x( x( op_index ) );
        } break;
        case Kinematics::InelasticElastic: {
          MX_ = mx_limits_.x( x( op_index ) );
          MY_ = event_->getOneByRole( Particle::IncomingBeam2 ).mass();
        } break;
        case Kinematics::InelasticInelastic: {
          MX_ = mx_limits_.x( x( op_index ) );
          MY_ = mx_limits_.x( x( op_index+1 ) );
        } break;
      }
      DebuggingInsideLoop( Form( "outgoing remnants invariant mass: %g / %g (%g < M(X/Y) < %g)", MX_, MY_, mx_limits_.min(), mx_limits_.max() ) );
    }

    void
    GenericKTProcess::computeIncomingFluxes( double x1, double q1t2, double x2, double q2t2 )
    {
      flux1_ = flux2_ = 0.;
      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic:
          flux1_ = elasticFlux( x1, q1t2 );
          flux2_ = elasticFlux( x2, q2t2 );
          break;
        case Kinematics::ElasticInelastic:
          flux1_ = elasticFlux( x1, q1t2 );
          flux2_ = inelasticFlux( x2, q2t2, MY_, cuts_.structure_functions );
          break;
        case Kinematics::InelasticElastic:
          flux1_ = inelasticFlux( x1, q1t2, MX_, cuts_.structure_functions );
          flux2_ = elasticFlux( x2, q2t2 );
          break;
        case Kinematics::InelasticInelastic:
          flux1_ = inelasticFlux( x1, q1t2, MX_, cuts_.structure_functions );
          flux2_ = inelasticFlux( x2, q2t2, MY_, cuts_.structure_functions );
          break;
        default: return;
      }
      flux1_ = std::max( flux1_, 1.e-20 );
      flux2_ = std::max( flux2_, 1.e-20 );
      DebuggingInsideLoop( Form( "Form factors: %g / %g", flux1_, flux2_ ) );
    }

    void
    GenericKTProcess::fillKinematics( bool )
    {
      fillPrimaryParticlesKinematics();
      fillCentralParticlesKinematics(); // process-dependent!
    }

    void
    GenericKTProcess::fillPrimaryParticlesKinematics()
    {
      //================================================================
      //     outgoing protons
      //================================================================
      Particle& op1 = event_->getOneByRole( Particle::OutgoingBeam1 ),
               &op2 = event_->getOneByRole( Particle::OutgoingBeam2 );

      op1.setMomentum( PX_ );
      op2.setMomentum( PY_ );

      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::ElasticInelastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::Unfragmented ); op2.setMass( MY_ );
          break;
        case Kinematics::InelasticElastic:
          op1.setStatus( Particle::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::InelasticInelastic:
          op1.setStatus( Particle::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::Unfragmented ); op2.setMass( MY_ );
          break;    
        default: { FatalError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" ); } break;
      }

      //================================================================
      //     incoming partons (photons, pomerons, ...)
      //================================================================
      //FIXME ensure the validity of this approach
      Particle& g1 = event_->getOneByRole( Particle::Parton1 ),
               &g2 = event_->getOneByRole( Particle::Parton2 );
      g1.setMomentum( event_->getOneByRole( Particle::IncomingBeam1 ).momentum()-PX_, true );
      g1.setStatus( Particle::Incoming );
      g2.setMomentum( event_->getOneByRole( Particle::IncomingBeam2 ).momentum()-PY_, true );
      g2.setStatus( Particle::Incoming );

      //================================================================
      //     two-parton system
      //================================================================
      event_->getOneByRole( Particle::Intermediate ).setMomentum( g1.momentum()+g2.momentum() );
    }

    double
    GenericKTProcess::elasticFlux( double x, double kt2 )
    {
      const double mp = ParticleProperties::mass( Proton ), mp2 = mp*mp;

      const double Q2_min = x*x*mp2/( 1.-x ), Q2_ela = ( kt2+x*x*mp2 )/( 1.-x );
      const FormFactors ff = FormFactors::ProtonElastic( Q2_ela );
      const double ela1 = ( 1.-x )*( 1.-Q2_min/Q2_ela );
      //const double ela3 = 1.-( Q2_ela-kt2 )/Q2_ela;

      const double f_ela = Constants::alphaEM/M_PI/kt2*( ela1*ff.FE + 0.5*x*x*ff.FM );

      // last factor below the Jacobian from dQ^2/Q^2 --> dkT^2/kT^2*(kT^2/Q^2)
      return f_ela*( 1.-x )*kt2 / Q2_ela;
    }

    double
    GenericKTProcess::inelasticFlux( double x, double kt2, double mx, const StructureFunctions::Type& sf )
    {
      const double mx2 = mx*mx, mp = ParticleProperties::mass( Proton ), mp2 = mp*mp;

      // F2 structure function
      const double Q2min = 1. / ( 1.-x )*( x*( mx2-mp2 ) + x*x*mp2 ),
                   Q2 = Q2min + kt2/( 1.-x );
      float xbj = Q2 / ( Q2 + mx2 - mp2 );

      const StructureFunctions str_fun = StructureFunctionsBuilder::get( sf, Q2, xbj );

      const double F1 = 0.5*( ( 1+4.*xbj*xbj*mp2/Q2 )*str_fun.F2 - str_fun.FL )/xbj;

      const double term1 = ( 1.-x )*( 1.-Q2min/Q2 );

      const double f_D = str_fun.F2/( mx2 + Q2 - mp2 ) * term1;
      const double f_C= 2.*F1/Q2;

      return Constants::alphaEM/M_PI*( 1.-x )/Q2*( f_D+0.5*x*x*f_C );
    }
  }
}
