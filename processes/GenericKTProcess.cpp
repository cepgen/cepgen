#include "GenericKTProcess.h"

GenericKTProcess::GenericKTProcess( const std::string& name_,
                                    const unsigned int& num_user_dimensions_,
                                    const Particle::ParticleCode& ip1_,
                                    const Particle::ParticleCode& op1_,
                                    const Particle::ParticleCode& ip2_,
                                    const Particle::ParticleCode& op2_) :
  GenericProcess( name_+" (kT-factorisation approach)" ),
  kNumUserDimensions( num_user_dimensions_ ),
  kIntermediatePart1( ip1_ ), kProducedPart1( op1_ )
{
  if ( ip2_==Particle::invalidParticle ) kIntermediatePart2 = kIntermediatePart1;
  if ( op2_==Particle::invalidParticle ) kProducedPart2 = kProducedPart1;
}

GenericKTProcess::~GenericKTProcess()
{}

void
GenericKTProcess::AddEventContent()
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
  GenericProcess::SetEventContent( is, os );
}

unsigned int
GenericKTProcess::GetNdim( const Kinematics::ProcessMode& process_mode_ ) const
{
  switch ( process_mode_ ) {
    default:
    case Kinematics::ElasticElastic:     return kNumRequiredDimensions+kNumUserDimensions;
    case Kinematics::ElasticInelastic:
    case Kinematics::InelasticElastic:   return kNumRequiredDimensions+kNumUserDimensions+1;
    case Kinematics::InelasticInelastic: return kNumRequiredDimensions+kNumUserDimensions+2;
  }
}

void
GenericKTProcess::AddPartonContent()
{
  // Incoming partons
  fQT1 = exp( fLogQmin+( fLogQmax-fLogQmin )*x( 0 ) );
  fQT2 = exp( fLogQmin+( fLogQmax-fLogQmin )*x( 1 ) );
  fPhiQT1 = 2.*Constants::Pi*x( 2 );
  fPhiQT2 = 2.*Constants::Pi*x( 3 );
  DebuggingInsideLoop( Form( "photons transverse virtualities (qt):\n\t"
                             "  mag = %f / %f (%.2f < log(qt) < %.2f)\n\t"
                             "  phi = %f / %f",
                             fQT1, fQT2, fLogQmin, fLogQmax, fPhiQT1, fPhiQT2 ) );
}

double
GenericKTProcess::ComputeWeight()
{
  AddPartonContent();
  PrepareKTKinematics();
  ComputeOutgoingPrimaryParticlesMasses();
  
  const double jac = ComputeJacobian(),
               integrand = ComputeKTFactorisedMatrixElement(),
               weight = jac*integrand;
  DebuggingInsideLoop( Form( "Jacobian = %f\n\tIntegrand = %f\n\tdW = %f", jac, integrand, weight ) );
  
  return weight;
}

void
GenericKTProcess::ComputeOutgoingPrimaryParticlesMasses()
{
  const unsigned int op_index = kNumRequiredDimensions+kNumUserDimensions;
  switch ( fCuts.kinematics ) {
    case Kinematics::ElectronProton: default: {
      InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
      exit( 0 ); } break;
    case Kinematics::ElasticElastic: 
      fMX = GetParticle( Particle::IncomingBeam1 )->M();
      fMY = GetParticle( Particle::IncomingBeam2 )->M();
      break;
    case Kinematics::ElasticInelastic:
      fMX = GetParticle( Particle::IncomingBeam1 )->M();
      fMY = fCuts.mxmin+( fCuts.mxmax-fCuts.mxmin )*x( op_index );
      break;
    case Kinematics::InelasticElastic:
      fMX = fCuts.mxmin+( fCuts.mxmax-fCuts.mxmin )*x( op_index );
      fMY = GetParticle( Particle::IncomingBeam2 )->M();
      break;
    case Kinematics::InelasticInelastic:
      fMX = fCuts.mxmin+( fCuts.mxmax-fCuts.mxmin )*x( op_index );
      fMY = fCuts.mxmin+( fCuts.mxmax-fCuts.mxmin )*x( op_index+1 );
      break;
  }
  DebuggingInsideLoop( Form( "outgoing remnants invariant mass: %f / %f (%.2f < M(X/Y) < %.2f)", fMX, fMY, fCuts.mxmin, fCuts.mxmax ) );
}

void
GenericKTProcess::ComputeIncomingFluxes( double x1, double q1t2, double x2, double q2t2 )
{
  fFlux1 = fFlux2 = 0.;
  switch ( fCuts.kinematics ) {
    case Kinematics::ElasticElastic:
      fFlux1 = PhotonFluxes::ProtonElastic( x1, q1t2 );
      fFlux2 = PhotonFluxes::ProtonElastic( x2, q2t2 );
      break;
    case Kinematics::ElasticInelastic:
      fFlux1 = PhotonFluxes::ProtonElastic( x1, q1t2 );
      fFlux2 = PhotonFluxes::ProtonInelastic( x2, q2t2, fMY );
      break;
    case Kinematics::InelasticElastic:
      fFlux1 = PhotonFluxes::ProtonInelastic( x1, q1t2, fMX );
      fFlux2 = PhotonFluxes::ProtonElastic( x2, q2t2 );
      break;
    case Kinematics::InelasticInelastic:
      fFlux1 = PhotonFluxes::ProtonInelastic( x1, q1t2, fMX );
      fFlux2 = PhotonFluxes::ProtonInelastic( x2, q2t2, fMY );
      break;
    default: return;
  }
  if ( fFlux1<1.e-20 ) fFlux1 = 0.;
  if ( fFlux2<1.e-20 ) fFlux2 = 0.;
  DebuggingInsideLoop( Form( "Form factors: %e / %e", fFlux1, fFlux2 ) );
}

void
GenericKTProcess::FillKinematics( bool )
{
  FillPrimaryParticlesKinematics();
  FillCentralParticlesKinematics();
}

void
GenericKTProcess::FillPrimaryParticlesKinematics()
{
  //=================================================================
  //     outgoing protons
  //=================================================================
  Particle *op1 = GetParticle( Particle::OutgoingBeam1 ),
           *op2 = GetParticle( Particle::OutgoingBeam2 );
  // off-shell particles (remnants?)
  bool os1 = false, os2 = false;
  switch ( fCuts.kinematics ) {
    case Kinematics::ElectronProton: default: {
      InError( "This kT factorisation process is intended for p-on-p collisions! Aborting!" );
      exit(0); } break;
    case Kinematics::ElasticElastic:
      op1->status = Particle::FinalState;
      op2->status = Particle::FinalState;
      break;
    case Kinematics::ElasticInelastic:
      op1->status = Particle::FinalState;
      op2->status = Particle::Undecayed; op2->SetM(); os2 = true;
      break;
    case Kinematics::InelasticElastic:
      op1->status = Particle::Undecayed; op1->SetM(); os1 = true;
      op2->status = Particle::FinalState;
      break;
    case Kinematics::InelasticInelastic:
      op1->status = Particle::Undecayed; op1->SetM(); os1 = true;
      op2->status = Particle::Undecayed; op2->SetM(); os2 = true;
      break;    
  }
  
  if ( !op1->SetMomentum( fPX, os1 ) ) { InError( Form( "Invalid outgoing proton 1: energy: %.2f", fPX.E() ) ); }
  if ( !op2->SetMomentum( fPY, os2 ) ) { InError( Form( "Invalid outgoing proton 2: energy: %.2f", fPY.E() ) ); }
  
  //=================================================================
  //     incoming partons (photons, pomerons, ...)
  //=================================================================
  //FIXME ensure the validity of this approach
  Particle *g1 = GetParticle( Particle::Parton1 ),
           *g2 = GetParticle( Particle::Parton2 );
  g1->SetMomentum( GetParticle( Particle::IncomingBeam1 )->GetMomentum()-fPX, true);
  g1->status = Particle::Incoming;
  g2->SetMomentum( GetParticle( Particle::IncomingBeam2 )->GetMomentum()-fPY, true);
  g2->status = Particle::Incoming;
}

double
GenericKTProcess::MinimalJacobian() const
{
  double jac = 1.;
  jac *= ( fLogQmax-fLogQmin )*fQT1; // d(q1t) . q1t
  jac *= ( fLogQmax-fLogQmin )*fQT2; // d(q2t) . q2t
  jac *= 2.*Constants::Pi; // d(phi1)
  jac *= 2.*Constants::Pi; // d(phi2)
  switch ( fCuts.kinematics ) {
    case Kinematics::ElasticElastic: default: break;
    case Kinematics::ElasticInelastic:   jac *= ( fCuts.mxmax-fCuts.mxmin )*2.*fMY; break;
    case Kinematics::InelasticElastic:   jac *= ( fCuts.mxmax-fCuts.mxmin )*2.*fMX; break;
    case Kinematics::InelasticInelastic: jac *= ( fCuts.mxmax-fCuts.mxmin )*2.*fMX;
                                         jac *= ( fCuts.mxmax-fCuts.mxmin )*2.*fMY; break;
  } // d(mx/y**2)
  return jac;
}
