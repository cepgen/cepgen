#include "GenericKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    GenericKTProcess::GenericKTProcess( const std::string& name_,
                                        const unsigned int& num_user_dimensions_,
                                        const Particle::ParticleCode& ip1,
                                        const Particle::ParticleCode& op1,
                                        const Particle::ParticleCode& ip2,
                                        const Particle::ParticleCode& op2) :
      GenericProcess( name_+" (kT-factorisation approach)" ),
      kNumUserDimensions( num_user_dimensions_ ),
      kIntermediatePart1( ip1 ), kProducedPart1( op1 )
    {
      if ( ip2==Particle::invalidParticle ) kIntermediatePart2 = kIntermediatePart1;
      if ( op2==Particle::invalidParticle ) kProducedPart2 = kProducedPart1;
    }

    GenericKTProcess::~GenericKTProcess()
    {}

    void
    GenericKTProcess::addEventContent()
    {
      IncomingState is; OutgoingState os;
      is.insert( ParticleWithRole( Particle::IncomingBeam1,    Particle::Proton ) );
      is.insert( ParticleWithRole( Particle::IncomingBeam2,    Particle::Proton ) );
      is.insert( ParticleWithRole( Particle::Parton1,          kIntermediatePart1 ) );
      is.insert( ParticleWithRole( Particle::Parton2,          kIntermediatePart2 ) );
      os.insert( ParticleWithRole( Particle::OutgoingBeam1,    Particle::Proton ) );
      os.insert( ParticleWithRole( Particle::OutgoingBeam2,    Particle::Proton ) );
      os.insert( ParticleWithRole( Particle::CentralParticle1, kProducedPart1 ) );
      os.insert( ParticleWithRole( Particle::CentralParticle2, kProducedPart2 ) );
      GenericProcess::setEventContent( is, os );
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
      switch ( cuts_.kinematics ) {
        case Kinematics::ElectronProton: default: {
          InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
          exit( 0 ); } break;
        case Kinematics::ElasticElastic: 
          MX_ = particlePtr( Particle::IncomingBeam1 )->mass();
          MY_ = particlePtr( Particle::IncomingBeam2 )->mass();
        break;
        case Kinematics::ElasticInelastic:
          MX_ = particlePtr( Particle::IncomingBeam1 )->mass();
          MY_ = cuts_.mx_min+( cuts_.mx_max-cuts_.mx_min )*x( op_index );
          break;
        case Kinematics::InelasticElastic:
          MX_ = cuts_.mx_min+( cuts_.mx_max-cuts_.mx_min )*x( op_index );
          MY_ = particlePtr( Particle::IncomingBeam2 )->mass();
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
      switch ( cuts_.kinematics ) {
        case Kinematics::ElasticElastic:
          flux1_ = PhotonFluxes::ProtonElastic( x1, q1t2 );
          flux2_ = PhotonFluxes::ProtonElastic( x2, q2t2 );
          break;
        case Kinematics::ElasticInelastic:
          flux1_ = PhotonFluxes::ProtonElastic( x1, q1t2 );
          flux2_ = PhotonFluxes::ProtonInelastic( x2, q2t2, MY_ );
          break;
        case Kinematics::InelasticElastic:
          flux1_ = PhotonFluxes::ProtonInelastic( x1, q1t2, MX_ );
          flux2_ = PhotonFluxes::ProtonElastic( x2, q2t2 );
          break;
        case Kinematics::InelasticInelastic:
          flux1_ = PhotonFluxes::ProtonInelastic( x1, q1t2, MX_ );
          flux2_ = PhotonFluxes::ProtonInelastic( x2, q2t2, MY_ );
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
      Particle *op1 = particlePtr( Particle::OutgoingBeam1 ),
               *op2 = particlePtr( Particle::OutgoingBeam2 );
      // off-shell particles (remnants?)
      bool os1 = false, os2 = false;
      switch ( cuts_.kinematics ) {
        case Kinematics::ElectronProton: default: {
          InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
          exit(0); } break;
        case Kinematics::ElasticElastic:
          op1->status = Particle::FinalState;
          op2->status = Particle::FinalState;
          break;
        case Kinematics::ElasticInelastic:
          op1->status = Particle::FinalState;
          op2->status = Particle::Undecayed; op2->setMass(); os2 = true;
          break;
        case Kinematics::InelasticElastic:
          op1->status = Particle::Undecayed; op1->setMass(); os1 = true;
          op2->status = Particle::FinalState;
          break;
        case Kinematics::InelasticInelastic:
          op1->status = Particle::Undecayed; op1->setMass(); os1 = true;
          op2->status = Particle::Undecayed; op2->setMass(); os2 = true;
          break;    
      }

      if ( !op1->setMomentum( PX_, os1 ) ) { InError( Form( "Invalid outgoing proton 1: energy: %.2f", PX_.energy() ) ); }
      if ( !op2->setMomentum( PY_, os2 ) ) { InError( Form( "Invalid outgoing proton 2: energy: %.2f", PY_.energy() ) ); }
  
      //=================================================================
      //     incoming partons (photons, pomerons, ...)
      //=================================================================
      //FIXME ensure the validity of this approach
      Particle *g1 = particlePtr( Particle::Parton1 ),
               *g2 = particlePtr( Particle::Parton2 );
      g1->setMomentum( particlePtr( Particle::IncomingBeam1 )->momentum()-PX_, true);
      g1->status = Particle::Incoming;
      g2->setMomentum( particlePtr( Particle::IncomingBeam2 )->momentum()-PY_, true);
      g2->status = Particle::Incoming;
    }

    double
    GenericKTProcess::minimalJacobian() const
    {
      double jac = 1.;
      jac *= ( log_qmax_-log_qmin_ )*qt1_; // d(q1t) . q1t
      jac *= ( log_qmax_-log_qmin_ )*qt2_; // d(q2t) . q2t
      jac *= 2.*M_PI; // d(phi1)
      jac *= 2.*M_PI; // d(phi2)
      switch ( cuts_.kinematics ) {
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
