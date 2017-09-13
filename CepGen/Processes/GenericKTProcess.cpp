#include "GenericKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    GenericKTProcess::GenericKTProcess( const std::string& name,
                                        const std::string& description,
                                        const unsigned int& num_user_dimensions,
                                        const std::array<Particle::ParticleCode,2>& partons,
                                        const std::vector<Particle::ParticleCode>& central ) :
      GenericProcess( name, description+" (kT-factorisation approach)" ),
      kNumUserDimensions( num_user_dimensions ),
      kIntermediateParts( partons ), kProducedParts( central )
    {}

    GenericKTProcess::~GenericKTProcess()
    {}

    void
    GenericKTProcess::addEventContent()
    {
      GenericProcess::setEventContent(
        { // incoming state
          { Particle::IncomingBeam1, Particle::Proton },
          { Particle::IncomingBeam2, Particle::Proton },
          { Particle::Parton1, kIntermediateParts[0] },
          { Particle::Parton2, kIntermediateParts[1] }
        },
        { // outgoing state
          { Particle::OutgoingBeam1, { Particle::Proton } },
          { Particle::OutgoingBeam2, { Particle::Proton } },
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
      qt1_ = exp( log_qmin_+( log_qmax_-log_qmin_ )*x( 0 ) );
      qt2_ = exp( log_qmin_+( log_qmax_-log_qmin_ )*x( 1 ) );
      phi_qt1_ = 2.*M_PI*x( 2 );
      phi_qt2_ = 2.*M_PI*x( 3 );
      DebuggingInsideLoop( Form( "photons transverse virtualities (qt):\n\t"
                                 "  mag = %f / %f (%.2f < log(qt) < %.2f)\n\t"
                                 "  phi = %f / %f",
                                 qt1_, qt2_, log_qmin_, log_qmax_, phi_qt1_, phi_qt2_ ) );
    }

    void
    GenericKTProcess::setKinematics( const Kinematics& kin )
    {
      cuts_ = kin;
      log_qmin_ = -10.; // FIXME //log_qmin_ = std::log( std::sqrt( cuts_.q2min ) );
      log_qmax_ = log( cuts_.initial_cuts[Cuts::qt].max() );
    }

    double
    GenericKTProcess::computeWeight()
    {
      addPartonContent();
      prepareKTKinematics();
      computeOutgoingPrimaryParticlesMasses();

      const double jac = computeJacobian(),
                   integrand = computeKTFactorisedMatrixElement(),
                   weight = jac*integrand;
      DebuggingInsideLoop( Form( "Jacobian = %f\n\tIntegrand = %f\n\tdW = %f", jac, integrand, weight ) );

      return weight;
    }

    void
    GenericKTProcess::computeOutgoingPrimaryParticlesMasses()
    {
      const unsigned int op_index = kNumRequiredDimensions+kNumUserDimensions;
      const Kinematics::Limits remn_mx_cuts = cuts_.remnant_cuts[Cuts::mass];
      switch ( cuts_.mode ) {
        case Kinematics::ElectronProton: default: {
          InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
          exit( 0 ); } break;
        case Kinematics::ElasticElastic: {
          MX_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
          MY_ = event_->getOneByRole( Particle::IncomingBeam2 ).mass();
        } break;
        case Kinematics::ElasticInelastic: {
          const double mx_min = remn_mx_cuts.min(), mx_range = remn_mx_cuts.range();
          MX_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
          MY_ = mx_min + mx_range*x( op_index );
        } break;
        case Kinematics::InelasticElastic: {
          const double mx_min = remn_mx_cuts.min(), mx_range = remn_mx_cuts.range();
          MX_ = mx_min + mx_range*x( op_index );
          MY_ = event_->getOneByRole( Particle::IncomingBeam2 ).mass();
        } break;
        case Kinematics::InelasticInelastic: {
          const double mx_min = remn_mx_cuts.min(), mx_range = remn_mx_cuts.range();
          MX_ = mx_min + mx_range*x( op_index );
          MY_ = mx_min + mx_range*x( op_index+1 );
        } break;
      }
      DebuggingInsideLoop( Form( "outgoing remnants invariant mass: %f / %f (%.2f < M(X/Y) < %.2f)", MX_, MY_, remn_mx_cuts.min(), remn_mx_cuts.max() ) );
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
          flux2_ = inelasticFlux( x2, q2t2, MY_ );
          break;
        case Kinematics::InelasticElastic:
          flux1_ = inelasticFlux( x1, q1t2, MX_ );
          flux2_ = elasticFlux( x2, q2t2 );
          break;
        case Kinematics::InelasticInelastic:
          flux1_ = inelasticFlux( x1, q1t2, MX_ );
          flux2_ = inelasticFlux( x2, q2t2, MY_ );
          break;
        default: return;
      }
      flux1_ = std::max( flux1_, 1.e-20 );
      flux2_ = std::max( flux2_, 1.e-20 );
      DebuggingInsideLoop( Form( "Form factors: %e / %e", flux1_, flux2_ ) );
    }

    void
    GenericKTProcess::fillKinematics( bool )
    {
      fillPrimaryParticlesKinematics();
      fillCentralParticlesKinematics();
    }

    void
    GenericKTProcess::fillPrimaryParticlesKinematics()
    {
      //=================================================================
      //     outgoing protons
      //=================================================================
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
          op2.setStatus( Particle::Undecayed ); op2.setMass( MY_ );
          break;
        case Kinematics::InelasticElastic:
          op1.setStatus( Particle::Undecayed ); op1.setMass( MX_ );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::InelasticInelastic:
          op1.setStatus( Particle::Undecayed ); op1.setMass( MX_ );
          op2.setStatus( Particle::Undecayed ); op2.setMass( MY_ );
          break;    
        default: { FatalError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" ); } break;
      }

      //=================================================================
      //     incoming partons (photons, pomerons, ...)
      //=================================================================
      //FIXME ensure the validity of this approach
      Particle& g1 = event_->getOneByRole( Particle::Parton1 ),
               &g2 = event_->getOneByRole( Particle::Parton2 );
      g1.setMomentum( event_->getOneByRole( Particle::IncomingBeam1 ).momentum()-PX_, true );
      g1.setStatus( Particle::Incoming );
      g2.setMomentum( event_->getOneByRole( Particle::IncomingBeam2 ).momentum()-PY_, true );
      g2.setStatus( Particle::Incoming );
    }

    double
    GenericKTProcess::minimalJacobian() const
    {
      double jac = 1.;
      jac *= ( log_qmax_-log_qmin_ )*qt1_; // d(q1t) . q1t
      jac *= ( log_qmax_-log_qmin_ )*qt2_; // d(q2t) . q2t
      jac *= 2.*M_PI; // d(phi1)
      jac *= 2.*M_PI; // d(phi2)

      const double mx_range = cuts_.remnant_cuts.at( Cuts::mass ).range();
      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic: default: break;
        case Kinematics::ElasticInelastic:   jac *= 2.* mx_range * MY_; break;
        case Kinematics::InelasticElastic:   jac *= 2.* mx_range * MX_; break;
        case Kinematics::InelasticInelastic: jac *= 2.* mx_range * MX_;
                                             jac *= 2.* mx_range * MY_; break;
      } // d(mx/y**2)
      return jac;
    }

    double
    GenericKTProcess::elasticFlux( double x, double kt2 )
    {
      const double mp = Particle::massFromPDGId( Particle::Proton ), mp2 = mp*mp;

      const double Q2_ela = ( kt2+x*x*mp2 )/( 1.-x );
      const FormFactors ela = FormFactors::ProtonElastic( Q2_ela );

      const double ela1 = pow( kt2/( kt2+x*x*mp2 ), 2 );
      const double ela2 = ela.FE;
      //const double ela3 = 1.-(Q2_ela-kt2)/Q2_ela;
      //const double ela3 = 1.-x*x*mp2/Q2_ela/(1.-x);
      //return Constants::alpha_em/M_PI*(1.-x+x*x/4.)*ela1*ela2*ela3/kt2;
      return Constants::alphaEM/M_PI*ela1*ela2/Q2_ela;
      //return Constants::alphaEM/M_PI*( ( 1.-x )*ela1*ela2*ela3 + x*x/2.*G_M*G_M )/kt2;
    }


    double
    GenericKTProcess::inelasticFlux( double x, double kt2, double mx )
    {
      const double mx2 = mx*mx, mp = Particle::massFromPDGId( Particle::Proton ), mp2 = mp*mp;

      // F2 structure function
      const double Q2min = 1. / ( 1.-x )*( x*( mx2-mp2 ) + x*x*mp2 ),
                   Q2 = Q2min + kt2 / ( 1.-x );
      float xbj = Q2 / ( Q2 + mx2 - mp2 );

      StructureFunctions sf;
      switch ( cuts_.structure_functions ) {
        case StructureFunctionsType::SzczurekUleshchenko:
          sf = SF::SzczurekUleshchenko( Q2, xbj ); break;
        case StructureFunctionsType::SuriYennie:
          sf = SF::SuriYennie( Q2, xbj ); break;
        case StructureFunctionsType::Fiore:
          sf = SF::FioreBrasse( Q2, xbj ); break;
        case StructureFunctionsType::ALLM91:
          sf = SF::ALLM( Q2, xbj, SF::ALLMParameterisation::allm91() ); break;
        case StructureFunctionsType::ALLM97:
          sf = SF::ALLM( Q2, xbj, SF::ALLMParameterisation::allm97() ); break;
        default: break; //FIXME
      }

      // Longitudinal/transverse virtual photon cross section R
      // from Sibirtsev and Blunden (Phys Rev C 88,065202 (2013))
      const double R = 0.014 * Q2 * ( exp( -0.07*Q2 ) + 41.*exp( -0.8*Q2 ) );
      const double F1 = sf.F2 * ( 1.+4*xbj*xbj*mp2/Q2 )/( 1.+R )/( 2.*xbj );

      const double term1 = ( 1.-x )*( 1.-Q2min/Q2 );
      const double f_D = sf.F2/( Q2 + mx2 - mp2 ) * term1, f_C= 2.*F1/Q2;

      return Constants::alphaEM/M_PI*( 1.-x )*kt2/Q2*( f_D + 0.5*x*x*f_C )/kt2;
    }
  }
}
