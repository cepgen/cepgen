#include "GenericKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    GenericKTProcess::GenericKTProcess( const std::string& name,
                                        const std::string& description,
                                        const unsigned int& num_user_dimensions,
                                        const Particle::ParticleCode& parton1,
                                        const Particle::ParticleCode& central1,
                                        const Particle::ParticleCode& parton2,
                                        const Particle::ParticleCode& central2 ) :
      GenericProcess( name, description+" (kT-factorisation approach)" ),
      kNumUserDimensions( num_user_dimensions ),
      kIntermediatePart1( parton1 ), kIntermediatePart2( parton2 ),
      kProducedPart1( central1 ), kProducedPart2( central2 )
    {
      if ( parton2  == Particle::invalidParticle ) kIntermediatePart2 = kIntermediatePart1;
      if ( central2 == Particle::invalidParticle ) kProducedPart2 = kProducedPart1;
    }

    GenericKTProcess::~GenericKTProcess()
    {}

    void
    GenericKTProcess::addEventContent()
    {
      GenericProcess::setEventContent(
        { // incoming state
          { Particle::IncomingBeam1, Particle::Proton },
          { Particle::IncomingBeam2, Particle::Proton },
          { Particle::Parton1, kIntermediatePart1 },
          { Particle::Parton2, kIntermediatePart2 }
        },
        { // outgoing state
          { Particle::OutgoingBeam1, { Particle::Proton } },
          { Particle::OutgoingBeam2, { Particle::Proton } },
          { Particle::CentralSystem, { kProducedPart1, kProducedPart2 } }
        }
      );
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
      switch ( cuts_.mode ) {
        case Kinematics::ElectronProton: default: {
          InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
          exit( 0 ); } break;
        case Kinematics::ElasticElastic: 
          MX_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
          MY_ = event_->getOneByRole( Particle::IncomingBeam2 ).mass();
        break;
        case Kinematics::ElasticInelastic:
          MX_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
          MY_ = cuts_.mx_min+( cuts_.mx_max-cuts_.mx_min )*x( op_index );
          break;
        case Kinematics::InelasticElastic:
          MX_ = cuts_.mx_min+( cuts_.mx_max-cuts_.mx_min )*x( op_index );
          MY_ = event_->getOneByRole( Particle::IncomingBeam2 ).mass();
          break;
        case Kinematics::InelasticInelastic:
          MX_ = cuts_.mx_min+( cuts_.mx_max-cuts_.mx_min )*x( op_index );
          MY_ = cuts_.mx_min+( cuts_.mx_max-cuts_.mx_min )*x( op_index+1 );
          break;
      }
      DebuggingInsideLoop( Form( "outgoing remnants invariant mass: %f / %f (%.2f < M(X/Y) < %.2f)", MX_, MY_, cuts_.mx_min, cuts_.mx_max ) );
    }

    void
    GenericKTProcess::computeIncomingFluxes( double x1, double q1t2, double x2, double q2t2 )
    {
      flux1_ = flux2_ = 0.;
      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic:
          flux1_ = Fluxes::Photon::ProtonElastic( x1, q1t2 );
          flux2_ = Fluxes::Photon::ProtonElastic( x2, q2t2 );
          break;
        case Kinematics::ElasticInelastic:
          flux1_ = Fluxes::Photon::ProtonElastic( x1, q1t2 );
          flux2_ = Fluxes::Photon::ProtonInelastic( x2, q2t2, MY_ );
          break;
        case Kinematics::InelasticElastic:
          flux1_ = Fluxes::Photon::ProtonInelastic( x1, q1t2, MX_ );
          flux2_ = Fluxes::Photon::ProtonElastic( x2, q2t2 );
          break;
        case Kinematics::InelasticInelastic:
          flux1_ = Fluxes::Photon::ProtonInelastic( x1, q1t2, MX_ );
          flux2_ = Fluxes::Photon::ProtonInelastic( x2, q2t2, MY_ );
          break;
        default: return;
      }
      if ( flux1_<1.e-20 ) flux1_ = 0.;
      if ( flux2_<1.e-20 ) flux2_ = 0.;
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
      // off-shell particles (remnants?)
      bool os1 = false, os2 = false;
      switch ( cuts_.mode ) {
        case Kinematics::ElectronProton: default: {
          InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
          exit(0); } break;
        case Kinematics::ElasticElastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::ElasticInelastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::Undecayed ); op2.setMass(); os2 = true;
          break;
        case Kinematics::InelasticElastic:
          op1.setStatus( Particle::Undecayed ); op1.setMass(); os1 = true;
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::InelasticInelastic:
          op1.setStatus( Particle::Undecayed ); op1.setMass(); os1 = true;
          op2.setStatus( Particle::Undecayed ); op2.setMass(); os2 = true;
          break;    
      }

      op1.setMomentum( PX_, os1 );
      op2.setMomentum( PY_, os2 );
  
      //=================================================================
      //     incoming partons (photons, pomerons, ...)
      //=================================================================
      //FIXME ensure the validity of this approach
      Particle& g1 = event_->getOneByRole( Particle::Parton1 ),
               &g2 = event_->getOneByRole( Particle::Parton2 );
      g1.setMomentum( event_->getOneByRole( Particle::IncomingBeam1 ).momentum()-PX_, true);
      g1.setStatus( Particle::Incoming );
      g2.setMomentum( event_->getOneByRole( Particle::IncomingBeam2 ).momentum()-PY_, true);
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
      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic: default: break;
        case Kinematics::ElasticInelastic:   jac *= ( cuts_.mx_max-cuts_.mx_min )*2.*MY_; break;
        case Kinematics::InelasticElastic:   jac *= ( cuts_.mx_max-cuts_.mx_min )*2.*MX_; break;
        case Kinematics::InelasticInelastic: jac *= ( cuts_.mx_max-cuts_.mx_min )*2.*MX_;
                                             jac *= ( cuts_.mx_max-cuts_.mx_min )*2.*MY_; break;
      } // d(mx/y**2)
      return jac;
    }
  }
}
