#include "GenericKTProcess.h"

GenericKTProcess::GenericKTProcess( const std::string& name_,
                                    const unsigned int& num_user_dimensions_,
                                    const Particle::ParticleCode& ip1_,
                                    const Particle::ParticleCode& op1_,
                                    const Particle::ParticleCode& ip2_,
                                    const Particle::ParticleCode& op2_) :
  GenericProcess( name_+" (kT-factorisation approach)" ),
  kNumRequiredDimensions( 4 ), kNumUserDimensions( num_user_dimensions_ ),
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

double
GenericKTProcess::ElasticFlux( double x_, double kt2_ ) const
{
  double f_ela;

  const double mp = Particle::GetMassFromPDGId(Particle::Proton),
               mp2 = mp*mp;

  const double Q2_ela = ( kt2_+x_*x_*mp2 )/( 1.-x_ );
  const FormFactors ela = ElasticFormFactors( Q2_ela, mp2 );
  
  /*const double G_dip = 1./pow(1.+Q2_ela/0.71, 2);
  const double G_E = G_dip;
  const double G_M = 2.79*G_dip;*/

  const double ela1 = pow( kt2_/( kt2_+x_*x_*mp2 ), 2 );
  const double ela2 = ela.FE;
  //const double ela3 = 1.-(Q2_ela-kt2_)/Q2_ela;
  //const double ela3 = 1.-pow(x_, 2)*mp2/Q2_ela/(1.-x_);
  //f_ela = alpha_em/Constants::Pi*(1.-x_+pow(x_, 2)/4.)*ela1*ela2*ela3/kt2_;
  f_ela = Constants::AlphaEM/Constants::Pi*ela1*ela2/Q2_ela;
  //f_ela = Constants::AlphaEM/Constants::Pi*((1.-x_)*ela1*ela2*ela3+pow(x_, 2)/2.*pow(G_M, 2))/kt2_;

  return f_ela;
}

#ifdef GRVPDF

double
GenericKTProcess::InelasticFlux( double x_, double kt2_, double mx_ ) const
{
  double f_ine;

  const double mx2 = mx_*mx_,
               mp = Particle::GetMassFromPDGId( Particle::Proton ),
               mp2 = mp*mp;
  //const double mpi = pow(Particle::GetMassFromPDGId(Particle::PiZero), 2);

  const double Q02 = 0.8; // introduced to shift the Q2 scale
  double term1, term2;
  double f_aux;

  // F2 structure function
  const double Q2min = 1. / ( 1.-x_ )*( x_*( mx2-mp2 ) + x_*x_*mp2 ),
               Q2 = kt2_ / ( 1.-x_ ) + Q2min;
  float x_Bjorken = Q2 / ( Q2+mx2-mp2 );

  float mu2 = Q2+Q02; // scale is shifted

  float xuv, xdv, xus, xds, xss, xg;
  grv95lo_( x_Bjorken, mu2, xuv, xdv, xus, xds, xss, xg );
  DebuggingInsideLoop( Form( "Form factor content at xB = %e (scale = %f GeV^2):\n\t"
                             "  valence quarks: u / d     = %e / %e\n\t"
                             "  sea quarks:     u / d / s = %e / %e / %e\n\t"
                             "  gluons:                   = %e",
                             x_Bjorken, mu2, xuv, xdv, xus, xds, xss, xg ) );

  const double F2_aux = 4./9.*( xuv + 2.*xus )
                      + 1./9.*( xdv + 2.*xds )
                      + 1./9.*2.*xss;

  /*F2_aux = 4./9.*(xuv + 2.*xus)
         + 1./9.*(0. + 2.*xds)
         + 1./9.*2.*xss;*/

  // F2 corrected for low Q^2 behaviour
  const double F2_corr = Q2 / ( Q2+Q02 ) * F2_aux;

  ///////term1 = pow(1.- x_/2.*(mx2-mp2+Q2)/Q2, 2);
  //term1 = (1.-x_*(mx2-mp2+Q2)/Q2);
  term1 = ( 1. - ( Q2-kt2_ ) / Q2 );
  //term1 = (1.-Q2min/Q2);
  //term1 = 1.;
  term2 = pow( kt2_ / ( kt2_+x_*(mx2-mp2)+x_*x_*mp2 ), 2 );

  f_aux = F2_corr/( mx2+Q2-mp2 )*term1*term2;

  f_ine = Constants::AlphaEM/Constants::Pi*( 1.-x_ )*f_aux/kt2_;

  return f_ine;
}

#endif
