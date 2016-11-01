#include "Particle.h"

Particle::Particle() :
  id( -1 ), charge( 1. ), name( "" ), role( UnknownRole ), status( Undefined ), helicity( 0. ),
  fMass( -1. ), fPDGid( invalidParticle ), fIsPrimary( true )
{}

Particle::Particle( Role role_, ParticleCode pdgId_ ) :
  id( -1 ), charge( 1. ), name( "" ), role( role_ ), status( Undefined ), helicity( 0. ),
  fMass( -1. ), fPDGid( pdgId_ ), fIsPrimary( true )
{
  if ( fPDGid!=invalidParticle ) {
    std::ostringstream o; o << fPDGid; name = o.str();
    SetM();
  }
}

Particle&
Particle::operator=( const Particle& part_ )
{
  fPDGid = part_.fPDGid;
  this->role = part_.role;
  if ( this->id==-1 ) this->id = part_.id;
  fMomentum = part_.fMomentum;
  this->SetM( part_.fMass );

  return *this;
}

bool
Particle::Valid()
{
  if ( fPDGid==invalidParticle ) return false;
  if ( fMomentum.P()==0. and M()==0. ) return false;
  return true;
}

std::string
Particle::GetLHEline( bool revert_ )
{
  std::stringstream line;

  if (revert_) fMomentum.SetP( 2, -fMomentum.P(2) );

  line << fPDGid << "\t";
  line << "1 1 2 0 0" << "\t";
  line << fMomentum.Px() << "\t";
  line << fMomentum.Py() << "\t";
  line << fMomentum.Pz() << "\t";
  line << fMomentum.E() << "\t";
  line << M() << "\t";
  line << "0." << "\t";
  line << "0."; //FIXME iz!!!
  return line.str();
}

bool
Particle::SetM( double m_ )
{
  double mass;
  if ( m_>=0. ) {
    fMass = m_;
    return true;
  }
  if ( fPDGid!=invalidParticle ) {
    mass = GetMassFromPDGId( fPDGid );
    if ( mass<0. ) return false;
    if ( fMomentum.E()<0. ) { // invalid energy
      fMass = mass;
      fMomentum.SetE( fMomentum.P2()+M2() );
      return true;
    }
    if ( E2()-fMomentum.P2()!=mass*mass ) {
      mass = std::sqrt( E2()-fMomentum.P2() );
    }
    fMass = mass;
    return true;
  }
  if ( E()>=0. and fMomentum.P()>=0. ) {
    fMass = sqrt( E2()-fMomentum.P2() );
    return true;
  }
  return false;
}

void
Particle::SetMother( Particle* part_ )
{
  fMothers.insert( part_->id );
  fIsPrimary = false;
  
  DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is the new mother of %2d (pdgId=%4d)",
                             part_->id+1, part_->GetPDGId(), id+1, fPDGid ) );
  
  part_->AddDaughter( this );
}

bool
Particle::AddDaughter( Particle* part_ )
{
  std::pair<ParticlesIds::iterator,bool> ret = fDaughters.insert( part_->id );

  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) {
    std::ostringstream os;
    for ( ParticlesIds::const_iterator it=fDaughters.begin(); it!=fDaughters.end(); it++) {
      os << Form("\n\t * id=%d", *it);
    }
    DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) has now %2d daughter(s):"
                               "%s", role, fPDGid, NumDaughters(), os.str().c_str() ) );
  }
  
  if ( ret.second ) {
    DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is a new daughter of %2d (pdgId=%4d)",
                               part_->role, part_->GetPDGId(), role, fPDGid ) );
    
    if ( !part_->Primary() and part_->GetMothersIds().size()<1 ) {
      part_->SetMother( this );
    }
  }

  return ret.second;
}

void
Particle::Dump() const
{  
  std::ostringstream osm, osd, os;
  if ( !Primary() ) {
    osm << ": mother(s): ";
    for ( ParticlesIds::const_iterator m=fMothers.begin(); m!=fMothers.end(); m++ ) {
      if ( m!=fMothers.begin() ) osm << ", ";
      osm << ( *m );
    }
  }
  const ParticlesIds daugh = GetDaughters();
  if ( daugh.size()!=0 ) {
    osd << ": id = ";
    for ( ParticlesIds::const_iterator it=daugh.begin(); it!=daugh.end(); it++ ) {
      if ( it!=daugh.begin() ) osd << ", ";
      osd << ( *it );
    }
  }
  os << " (" << fPDGid << ")";
  if ( os.str()==" ()" ) os.str("");
  Information( Form(
    "Dumping a particle with id=%3d, role=%3d, status=% 3d\n\t"
    "PDG Id:%4d%s, mass = %5.4f GeV\n\t"
    "(E,P) = (%4.2f, %4.2f, %4.2f, %4.2f) GeV\t"
    "(|P| = p = %4.2f GeV)\n\t"
    " Pt = %5.4f GeV, eta = %4.3f, phi = % 4.3f\n\t"
    "Primary? %s%s\n\t"
    "%d daughter(s)%s",
    id, role, status, fPDGid, os.str().c_str(),
    M(), E(), fMomentum.Px(), fMomentum.Py(), fMomentum.Pz(),
    fMomentum.P(), fMomentum.Pt(), fMomentum.Eta(), fMomentum.Phi(),
    yesno( Primary() ), osm.str().c_str(), NumDaughters(), osd.str().c_str() )
  );
}

void
Particle::LorentzBoost( double m_, const Particle::Momentum& mom_ )
{
  double pf4, fn;

  if (mom_.P(3)!=m_) {
    pf4 = 0.;
    for ( unsigned int i=0; i<4; i++ ) {
      pf4 += fMomentum.P(i)*mom_.P(i);
    }
    pf4 /= m_;
    fn = ( pf4+E() )/( fMomentum.P(3)+m_ );
    for ( unsigned int i=0; i<3; i++ ) {
      fMomentum.SetP( i, fMomentum.P(i)+fn*mom_.P(i) );
    }
  }
}

double*
Particle::LorentzBoost( const Particle::Momentum& mom_ )
{
  double p2, gamma, bp, gamma2;

  p2 = mom_.P2();
  gamma = 1./sqrt( 1.-p2 );
  bp = 0.;
  for ( unsigned int i=0; i<3; i++ ) bp+= mom_.P(i)*fMomentum.P(i);

  if ( p2>0. ) gamma2 = (gamma-1.)/p2;
  else gamma2 = 0.;

  for ( unsigned int i=0; i<3; i++ ) {
    __tmp3[i] = fMomentum.P(i) + gamma2*bp*mom_.P(i)+gamma*mom_.P(i)*E();
  }
  return __tmp3;
}

double
Particle::GetMassFromPDGId( const Particle::ParticleCode& pdgId_ )
{
  switch (pdgId_) {
    case dQuark:       return 0.33;           // mass from PYTHIA6.4
    case uQuark:       return 0.33;           // mass from PYTHIA6.4
    case Electron:     return 0.510998928e-3;
    case ElectronNeutrino: return 0.;
    case Muon:         return 0.1056583715;
    case MuonNeutrino: return 0.;
    case Tau:          return 1.77682;
    case TauNeutrino:  return 0.;
    case Gluon:        return 0.;
    case Z:            return 91.1876;
    case WPlus:        return 80.385;
    case Photon:       return 0.;
    case PiPlus:       return 0.13957018;
    case PiZero:       return 0.1349766;
    case JPsi:         return 20.;            //FIXME FIXME FIXME
    case ud0Diquark:   return 0.57933;
    case ud1Diquark:   return 0.77133;
    case uu1Diquark:   return 0.77133;
    case Proton:       return 0.938272046;
    case Neutron:      return 0.939565346;
    case Upsilon1S:    return 9.46030;
    case Upsilon2S:    return 10.02326;
    case Upsilon3S:    return 10.3552;
    case Rho770_0:     return 0.77526;
    case Rho1450_0:    return 1.465;
    case Rho1700_0:    return 1.720;
    case h1380_1:      return 1.38619;
    case invalidParticle:
    default:           return -1.;
  }
}

double
Particle::GetWidthFromPDGId( const Particle::ParticleCode& pdgId_ )
{
  switch (pdgId_) {
    case JPsi:      return 5.; //FIXME
    case Z:         return 2.4952;
    case WPlus:     return 2.085;
    case Upsilon1S: return 54.02e-6;
    case Upsilon2S: return 31.98e-6;
    case Upsilon3S: return 20.32e-6;
    case Rho770_0:  return 0.150; // PDG
    case Rho1450_0: return 0.400; // PDG
    case Rho1700_0: return 0.250; // PDG
    default:        return -1.;
  }
}

std::ostream&
operator<<( std::ostream& os, const Particle::ParticleCode& pc )
{
  switch (pc) {
    case Particle::dQuark:       os << "d quark"; break;
    case Particle::uQuark:       os << "u quark"; break;
    case Particle::Electron:     os << "electron"; break;
    case Particle::ElectronNeutrino: os << "electron neutrino"; break;
    case Particle::Muon:         os << "muon"; break;
    case Particle::MuonNeutrino: os << "muon neutrino"; break;
    case Particle::Tau:          os << "tau"; break;
    case Particle::TauNeutrino:  os << "tau neutrino"; break;
    case Particle::Gluon:        os << "gluon"; break;
    case Particle::Photon:       os << "photon"; break;
    case Particle::Z:            os << "Z"; break;
    case Particle::WPlus:        os << "W+"; break;
    case Particle::PiPlus:       os << "pi+"; break;
    case Particle::PiZero:       os << "pi0"; break;
    case Particle::Rho770_0:     os << "rho(770)0"; break;
    case Particle::Rho1450_0:    os << "rho(1450)0"; break;
    case Particle::Rho1700_0:    os << "rho(1700)0"; break;
    case Particle::h1380_1:      os << "h(1380)1"; break;
    case Particle::Omega782:     os << "omega(782)"; break;
    case Particle::JPsi:         os << "J/Psi"; break;
    case Particle::Phi1680:      os << "phi(1680)"; break;
    case Particle::Upsilon1S:    os << "Upsilon(1S)"; break;
    case Particle::Upsilon2S:    os << "Upsilon(2S)"; break;
    case Particle::Upsilon3S:    os << "Upsilon(3S)"; break;
    case Particle::ud0Diquark:   os << "(ud)0 di-quark"; break;
    case Particle::ud1Diquark:   os << "(ud)1 di-quark"; break;
    case Particle::uu1Diquark:   os << "(uu)1 di-quark"; break;
    case Particle::Proton:       os << "proton"; break;
    case Particle::Neutron:      os << "neutron"; break;
    case Particle::Pomeron:      os << "pomeron"; break;
    case Particle::Reggeon:      os << "reggeon"; break;
    case Particle::invalidParticle: os << "<invalid>"; break;
  }
  return os;
}

//----- Particle momentum methods

Particle::Momentum&
Particle::Momentum::operator+=( const Particle::Momentum& mom_ )
{
  fPx += mom_.fPx;
  fPy += mom_.fPy;
  fPz += mom_.fPz;
  fE += mom_.fE; //FIXME not supposed to be this way!
  ComputeP();
  return *this;
}

Particle::Momentum&
Particle::Momentum::operator-=( const Particle::Momentum& mom_ )
{
  fPx -= mom_.fPx;
  fPy -= mom_.fPy;
  fPz -= mom_.fPz;
  fE -= mom_.fE; //FIXME not supposed to be this way!
  ComputeP();
  return *this;
}

void
Particle::Momentum::operator=( const Particle::Momentum& mom_ )
{
  fPx = mom_.fPx; fPy = mom_.fPy; fPz = mom_.fPz; fP = mom_.fP;
  fE = mom_.fE;
}

double
Particle::Momentum::ThreeProduct( const Particle::Momentum& mom_ ) const
{
  DebuggingInsideLoop( Form( "  (%f, %f, %f, %f)\n\t* (%f, %f, %f, %f)\n\t= %f",
    fPx, fPy, fPz, fE,
    mom_.fPx, mom_.fPy, mom_.fPz, mom_.fE,
    fPx*mom_.fPx+fPy*mom_.fPy+fPz*mom_.fPz
  ) );
  return fPx*mom_.fPx+fPy*mom_.fPy+fPz*mom_.fPz;
}

double
Particle::Momentum::FourProduct( const Particle::Momentum& mom_ ) const
{
  DebuggingInsideLoop( Form( "  (%f, %f, %f, %f)\n\t* (%f, %f, %f, %f)\n\t= %f",
    fPx, fPy, fPz, fE,
    mom_.fPx, mom_.fPy, mom_.fPz, mom_.fE,
    fPx*mom_.fPx+fPy*mom_.fPy+fPz*mom_.fPz
  ) );
  return fE*mom_.fE-ThreeProduct(mom_);
}

double
Particle::Momentum::operator*=( const Particle::Momentum& mom_ )
{
  return ThreeProduct( mom_ );
}

Particle::Momentum&
Particle::Momentum::operator*=( double c )
{
  fPx *= c;
  fPy *= c;
  fPz *= c;
  ComputeP();
  return *this;
}

Particle::Momentum
operator*( const Particle::Momentum& mom, double c )
{
  Particle::Momentum out = mom;
  out *= c;
  return out;
}

Particle::Momentum
operator*( double c, const Particle::Momentum& mom )
{
  Particle::Momentum out = mom;
  out *= c;
  return out;
}

Particle::Momentum
operator+( const Particle::Momentum& mom1, const Particle::Momentum& mom2 )
{
  Particle::Momentum out = mom1;
  out += mom2;
  return out;
}

Particle::Momentum
operator-( const Particle::Momentum& mom1, const Particle::Momentum& mom2 )
{
  Particle::Momentum out = mom1;
  out -= mom2;
  return out;
}

double
operator*( const Particle::Momentum& mom1, const Particle::Momentum& mom2 )
{
  Particle::Momentum tmp = mom1;
  return tmp.ThreeProduct( mom2 );
}

void
Particle::Momentum::BetaGammaBoost( double gamma, double betagamma )
{
  const double pz = fPz, e = fE;
  fPz = gamma*pz+betagamma*e;
  fE  = gamma*e +betagamma*pz;
  ComputeP();
}

void
Particle::Momentum::LorentzBoost( const Particle::Momentum& p )
{
  const double m = p.M();
  if ( m==p.E() ) return;
  
  Particle::Momentum mom_old = *this;
  const double pf4 = mom_old*p/m,
               fn = ( pf4+mom_old.E() )/( p.E()+m );
  Particle::Momentum mom_new = mom_old-p*fn; mom_new.SetE(pf4);
  SetMomentum( mom_new );
}

void
Particle::Momentum::RotatePhi( double phi, double sign )
{
  const double px = fPx*cos( phi )+fPy*sin( phi )*sign,
               py =-fPx*sin( phi )+fPy*cos( phi )*sign;
  fPx = px;
  fPy = py;
}

void
Particle::Momentum::RotateThetaPhi( double theta_, double phi_ )
{
  double rotmtx[3][3], mom[3]; //FIXME check this! cos(phi)->-sin(phi) & sin(phi)->cos(phi) --> phi->phi+pi/2 ?
  rotmtx[0][0] = -sin( phi_ ); rotmtx[0][1] = -cos( theta_ )*cos( phi_ ); rotmtx[0][2] =  sin( theta_ )*cos( phi_ );
  rotmtx[1][0] =  cos( phi_ ); rotmtx[1][1] = -cos( theta_ )*sin( phi_ ); rotmtx[1][2] =  sin( theta_ )*sin( phi_ );
  rotmtx[2][0] =  0.;          rotmtx[2][1] =  sin( theta_ );             rotmtx[2][2] =  cos( theta_ );

  for (int i=0; i<3; i++) {
    mom[i] = 0.;
    for (int j=0; j<3; j++) {
      mom[i] += rotmtx[i][j]*P( j );
    }
  }

  SetP( mom[0], mom[1], mom[2] );
}

std::ostream&
operator<<( std::ostream& os, const Particle::Momentum& mom )
{
  os << "(E, p) = (" << mom.fE << ", " << mom.fPx << ", " << mom.fPy << ", " << mom.fPz << ")";
  return os;
}

