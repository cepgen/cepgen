#include "GenericProcess.h"

GenericProcess::GenericProcess( const std::string& name_ ) :
  fX( 0 ), fNumDimensions( 0 ), fEvent( new Event ),
  fIsPointSet( false ), fIsInStateSet( false ), fIsOutStateSet( false ), fIsKinematicSet( false ),
  fName( name_ ),
  fTotalGenTime( 0. ), fNumGenEvents( 0 )
{}

GenericProcess::~GenericProcess()
{
  if ( fIsPointSet ) delete[] fX;
  delete fEvent;
}

void
GenericProcess::SetPoint( const unsigned int ndim_, double x_[] )
{
  // Number of dimensions on which the integration will be performed
  fNumDimensions = ndim_;

  // Phase space coordinate becomes a protected attribute
  if ( !fX ) fX = new double[ndim_];

  std::copy( x_, x_+ndim_, fX );  
  fIsPointSet = true;
  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) { DumpPoint( DebugMessage ); }
}

void
GenericProcess::PrepareKinematics()
{
  if ( !IsKinematicsDefined() ) return; // FIXME dump some information...
  fSqS = CMEnergy( *GetParticle( Particle::IncomingBeam1 ),
                   *GetParticle( Particle::IncomingBeam2 ) );
  fS = fSqS*fSqS;
  
  fW1 = GetParticle( Particle::IncomingBeam1 )->M2();
  fW2 = GetParticle( Particle::IncomingBeam2 )->M2();

  Debugging( Form( "Kinematics successfully prepared! sqrt(s) = %.2f", fSqS ) );
}

void
GenericProcess::DumpPoint( const ExceptionType& et=Information )
{
  std::ostringstream os;
  for ( unsigned int i=0; i<fNumDimensions; i++ ) {
    os << Form( "  x(%2d) = %8.6f\n\t", i, fX[i] );
  }
  if ( et<DebugMessage ) { Information( Form( "Number of integration parameters: %d\n\t"
                                              "%s", fNumDimensions, os.str().c_str() ) ); }
  else                   { Debugging( Form( "Number of integration parameters: %d\n\t"
                                            "%s", fNumDimensions, os.str().c_str() ) ); }
}

void
GenericProcess::SetEventContent( const IncomingState& is, const OutgoingState& os )
{  
  for ( IncomingState::const_iterator ip=is.begin(); ip!=is.end(); ip++ ) { fEvent->AddParticle( Particle( ip->first, ip->second ) ); }

  // Prepare the central system if not already there
  IncomingState::const_iterator central_system = is.find( Particle::CentralSystem );
  if ( central_system==is.end() ) {
    Particle* moth = GetParticle( Particle::Parton1 );
    Particle cs( Particle::CentralSystem, moth->GetPDGId() );
    cs.SetMother( moth );
    fEvent->AddParticle( cs );
  }

  for ( OutgoingState::const_iterator op=os.begin(); op!=os.end(); op++ ) { fEvent->AddParticle( Particle( op->first, op->second ) ); }
  
  // Incoming particles (incl. eventual partons)
  for ( IncomingState::const_iterator ip=is.begin(); ip!=is.end(); ip++ ) {
    Particle* p = GetParticle( ip->first );
    p->status = Particle::Undefined;
    switch ( ip->first ) {
      case Particle::IncomingBeam1:
      case Particle::IncomingBeam2: break;
      case Particle::Parton1:       p->SetMother( GetParticle( Particle::IncomingBeam1 ) ); break;
      case Particle::Parton2:       p->SetMother( GetParticle( Particle::IncomingBeam2 ) ); break;
      case Particle::CentralSystem: p->SetMother( GetParticle( Particle::Parton1 ) ); break;
      default: break;
    }
  }
  // Outgoing particles (central, and outgoing primary particles or remnants)
  for ( OutgoingState::const_iterator op=os.begin(); op!=os.end(); op++ ) {
    Particle* p = GetParticle( op->first );
    p->status = Particle::Undefined;
    switch ( op->first ) {
      case Particle::OutgoingBeam1:    p->SetMother( GetParticle( Particle::IncomingBeam1 ) ); break;
      case Particle::OutgoingBeam2:    p->SetMother( GetParticle( Particle::IncomingBeam2 ) ); break;
      case Particle::CentralParticle1: p->SetMother( GetParticle( Particle::CentralSystem ) ); break;
      case Particle::CentralParticle2: p->SetMother( GetParticle( Particle::CentralSystem ) ); break;
      default: break;
    }
  }
  fEvent->Init();
}

void
GenericProcess::SetIncomingKinematics( const Particle::Momentum& p1, const Particle::Momentum& p2 )
{
  if ( !GetParticle( Particle::IncomingBeam1 )->SetMomentum( p1 ) ) { InError( "Invalid incoming beam 1" ); }
  if ( !GetParticle( Particle::IncomingBeam2 )->SetMomentum( p2 ) ) { InError( "Invalid incoming beam 2" ); }
}

void
GenericProcess::GetFormFactors( double q1, double q2, FormFactors& fp1, FormFactors& fp2 ) const
{
  const double mx2 = fMX*fMX, my2 = fMY*fMY;

  bool inel_p1 = false,
       inel_p2 = false;

  switch ( fCuts.kinematics ) {
    case Kinematics::ElectronElectron: {
      fp1 = TrivialFormFactors(); // electron (trivial) form factor
      fp2 = TrivialFormFactors(); // electron (trivial) form factor
    } break;
    case Kinematics::ProtonElectron: {
      fp1 = ElasticFormFactors( -fT1, fW1 ); // proton elastic form factor
      fp2 = TrivialFormFactors(); // electron (trivial) form factor
    } break;
    case Kinematics::ElectronProton: {
      fp1 = TrivialFormFactors(); // electron (trivial) form factor
      fp2 = ElasticFormFactors( -fT2, fW2 ); // proton elastic form factor
    } break;
    case Kinematics::ElasticElastic: {
      fp1 = ElasticFormFactors( -fT1, fW1 ); // proton elastic form factor
      fp2 = ElasticFormFactors( -fT2, fW2 ); // proton elastic form factor
    } break;
    case Kinematics::ElasticInelastic: {
      fp1 = ElasticFormFactors( -fT1, fW1 );
      inel_p2 = true;
    } break;
    case Kinematics::InelasticElastic: {
      inel_p1 = true;
      fp2 = ElasticFormFactors( -fT2, fW2 );
    } break;
    case Kinematics::InelasticInelastic: {
      inel_p1 = inel_p2 = true;
    } break;
  }
  switch ( fCuts.remnant_mode ) {
    case SuriYennie:
    default: {
      if ( inel_p1 ) fp1 = SuriYennieFormFactors( -fT1, fW1, mx2 );
      if ( inel_p2 ) fp2 = SuriYennieFormFactors( -fT2, fW2, my2 );
    } break;
    case Fiore:
    case FioreSea:
    case FioreVal: { // low-Q2 inelastic form factor
      if ( inel_p1 ) fp1 = FioreBrasseFormFactors( -fT1, fW1, mx2 );
      if ( inel_p2 ) fp2 = FioreBrasseFormFactors( -fT2, fW2, my2 );
    } break;
    case SzczurekUleshchenko: {
      if ( inel_p1 ) fp1 = SzczurekUleshchenkoFormFactors( -fT1, fW1, mx2 );
      if ( inel_p2 ) fp2 = SzczurekUleshchenkoFormFactors( -fT2, fW2, my2 );
    } break;
  }
}

std::ostream&
operator<<( std::ostream& os, const GenericProcess& proc )
{
  os << proc.GetName().c_str();
  return os;
}

std::ostream&
operator<<( std::ostream& os, const GenericProcess* proc )
{
  os << proc->GetName().c_str();
  return os;
}
