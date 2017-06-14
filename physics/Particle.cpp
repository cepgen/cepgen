#include "Particle.h"

Particle::Particle() :
  id( -1 ), charge( 1. ), name( "" ), role( UnknownRole ), status( Undefined ), helicity( 0. ),
  mass_( -1. ), pdg_id_( invalidParticle ), is_primary_( true )
{}

Particle::Particle( Role role, ParticleCode pdgId ) :
  id( -1 ), charge( 1. ), name( "" ), role( role ), status( Undefined ), helicity( 0. ),
  mass_( -1. ), pdg_id_( pdgId ), is_primary_( true )
{
  if ( pdg_id_!=invalidParticle ) {
    std::ostringstream o; o << pdg_id_; name = o.str();
    setMass();
  }
}

Particle&
Particle::operator=( const Particle& part )
{
  pdg_id_ = part.pdg_id_;
  this->role = part.role;
  if ( this->id==-1 ) this->id = part.id;
  momentum_ = part.momentum_;
  this->setMass( part.mass_ );

  return *this;
}

bool
Particle::valid()
{
  if ( pdg_id_==invalidParticle ) return false;
  if ( momentum_.p()==0. and mass()==0. ) return false;
  return true;
}

std::string
Particle::lheLine( bool revert )
{
  std::stringstream line;

  if ( revert ) momentum_.setP( 2, -momentum_.p( 2 ) );

  line << pdg_id_ << "\t";
  line << "1 1 2 0 0" << "\t";
  line << momentum_.px() << "\t";
  line << momentum_.py() << "\t";
  line << momentum_.pz() << "\t";
  line << momentum_.energy() << "\t";
  line << mass() << "\t";
  line << "0." << "\t";
  line << "0."; //FIXME iz!!!
  return line.str();
}

bool
Particle::setMass( double m )
{
  double mass;
  if ( m>=0. ) {
    mass_ = m;
    return true;
  }
  if ( pdg_id_!=invalidParticle ) {
    mass = massFromPDGId( pdg_id_ );
    if ( mass<0. ) return false;
    if ( momentum_.energy()<0. ) { // invalid energy
      mass_ = mass;
      momentum_.setEnergy( momentum_.p2()+mass2() );
      return true;
    }
    if ( energy2()-momentum_.p2()!=mass*mass ) {
      mass = std::sqrt( energy2()-momentum_.p2() );
    }
    mass_ = mass;
    return true;
  }
  if ( energy()>=0. and momentum_.p()>=0. ) {
    mass_ = sqrt( energy2()-momentum_.p2() );
    return true;
  }
  return false;
}

void
Particle::setMother( Particle* part )
{
  mothers_.insert( part->id );
  is_primary_ = false;
  
  DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is the new mother of %2d (pdgId=%4d)",
                             part->id+1, part->pdgId(), id+1, pdg_id_ ) );
  
  part->addDaughter( this );
}

bool
Particle::addDaughter( Particle* part )
{
  std::pair<ParticlesIds::iterator,bool> ret = daughters_.insert( part->id );

  if ( Logger::GetInstance()->Level>=Logger::DebugInsideLoop ) {
    std::ostringstream os;
    for ( ParticlesIds::const_iterator it=daughters_.begin(); it!=daughters_.end(); it++) {
      os << Form("\n\t * id=%d", *it);
    }
    DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) has now %2d daughter(s):"
                               "%s", role, pdg_id_, numDaughters(), os.str().c_str() ) );
  }
  
  if ( ret.second ) {
    DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is a new daughter of %2d (pdgId=%4d)",
                               part->role, part->pdgId(), role, pdg_id_ ) );
    
    if ( !part->primary() and part->mothersIds().size()<1 ) {
      part->setMother( this );
    }
  }

  return ret.second;
}

void
Particle::dump() const
{  
  std::ostringstream osm, osd, os;
  if ( !primary() ) {
    osm << ": mother(s): ";
    for ( ParticlesIds::const_iterator m=mothers_.begin(); m!=mothers_.end(); m++ ) {
      if ( m!=mothers_.begin() ) osm << ", ";
      osm << ( *m );
    }
  }
  const ParticlesIds daugh = daughters();
  if ( daugh.size()!=0 ) {
    osd << ": id = ";
    for ( ParticlesIds::const_iterator it=daugh.begin(); it!=daugh.end(); it++ ) {
      if ( it!=daugh.begin() ) osd << ", ";
      osd << ( *it );
    }
  }
  os << " (" << pdg_id_ << ")";
  if ( os.str()==" ()" ) os.str("");
  Information( Form(
    "Dumping a particle with id=%3d, role=%3d, status=% 3d\n\t"
    "PDG Id:%4d%s, mass = %5.4f GeV\n\t"
    "(E,P) = (%4.2f, %4.2f, %4.2f, %4.2f) GeV\t"
    "(|P| = p = %4.2f GeV)\n\t"
    " Pt = %5.4f GeV, eta = %4.3f, phi = % 4.3f\n\t"
    "Primary? %s%s\n\t"
    "%d daughter(s)%s",
    id, role, status, pdg_id_, os.str().c_str(),
    mass(), energy(), momentum_.px(), momentum_.py(), momentum_.pz(),
    momentum_.p(), momentum_.pt(), momentum_.eta(), momentum_.phi(),
    yesno( primary() ), osm.str().c_str(), numDaughters(), osd.str().c_str() )
  );
}

void
Particle::lorentzBoost( double m, const Particle::Momentum& mom )
{
  double pf4, fn;

  if ( mom.p( 3 )!=m ) {
    pf4 = 0.;
    for ( unsigned int i=0; i<4; i++ ) {
      pf4 += momentum_.p( i )*mom.p( i );
    }
    pf4 /= m;
    fn = ( pf4+energy() )/( momentum_.p( 3 )+m );
    for ( unsigned int i=0; i<3; i++ ) {
      momentum_.setP( i, momentum_.p( i )+fn*mom.p( i ) );
    }
  }
}

double*
Particle::lorentzBoost( const Particle::Momentum& mom )
{
  double p2, gamma, bp, gamma2;

  p2 = mom.p2();
  gamma = 1./sqrt( 1.-p2 );
  bp = 0.;
  for ( unsigned int i=0; i<3; i++ ) bp+= mom.p(i)*momentum_.p(i);

  if ( p2>0. ) gamma2 = (gamma-1.)/p2;
  else gamma2 = 0.;

  for ( unsigned int i=0; i<3; i++ ) {
    __tmp3[i] = momentum_.p(i) + gamma2*bp*mom.p(i)+gamma*mom.p(i)*energy();
  }
  return __tmp3;
}

double
Particle::massFromPDGId( const Particle::ParticleCode& pdgId_ )
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
Particle::widthFromPDGId( const Particle::ParticleCode& pdgId_ )
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
Particle::Momentum::operator+=( const Particle::Momentum& mom )
{
  px_ += mom.px_;
  py_ += mom.py_;
  pz_ += mom.pz_;
  energy_ += mom.energy_; //FIXME not supposed to be this way!
  computeP();
  return *this;
}

Particle::Momentum&
Particle::Momentum::operator-=( const Particle::Momentum& mom )
{
  px_ -= mom.px_;
  py_ -= mom.py_;
  pz_ -= mom.pz_;
  energy_ -= mom.energy_; //FIXME not supposed to be this way!
  computeP();
  return *this;
}

void
Particle::Momentum::operator=( const Particle::Momentum& mom )
{
  px_ = mom.px_; py_ = mom.py_; pz_ = mom.pz_; p_ = mom.p_;
  energy_ = mom.energy_;
}

double
Particle::Momentum::threeProduct( const Particle::Momentum& mom ) const
{
  DebuggingInsideLoop( Form( "  (%f, %f, %f, %f)\n\t* (%f, %f, %f, %f)\n\t= %f",
    px_, py_, pz_, energy_,
    mom.px_, mom.py_, mom.pz_, mom.energy_,
    px_*mom.px_+py_*mom.py_+pz_*mom.pz_
  ) );
  return px_*mom.px_+py_*mom.py_+pz_*mom.pz_;
}

double
Particle::Momentum::fourProduct( const Particle::Momentum& mom ) const
{
  DebuggingInsideLoop( Form( "  (%f, %f, %f, %f)\n\t* (%f, %f, %f, %f)\n\t= %f",
    px_, py_, pz_, energy_,
    mom.px_, mom.py_, mom.pz_, mom.energy_,
    px_*mom.px_+py_*mom.py_+pz_*mom.pz_
  ) );
  return energy_*mom.energy_-threeProduct(mom);
}

double
Particle::Momentum::operator*=( const Particle::Momentum& mom )
{
  return threeProduct( mom );
}

Particle::Momentum&
Particle::Momentum::operator*=( double c )
{
  px_ *= c;
  py_ *= c;
  pz_ *= c;
  computeP();
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
  return tmp.threeProduct( mom2 );
}

void
Particle::Momentum::betaGammaBoost( double gamma, double betagamma )
{
  const double pz = pz_, e = energy_;
  pz_ = gamma*pz+betagamma*e;
  energy_  = gamma*e +betagamma*pz;
  computeP();
}

void
Particle::Momentum::lorentzBoost( const Particle::Momentum& p )
{
  const double m = p.mass();
  if ( m==p.energy() ) return;
  
  Particle::Momentum mom_old = *this;
  const double pf4 = mom_old*p/m,
               fn = ( pf4+mom_old.energy() )/( p.energy()+m );
  Particle::Momentum mom_new = mom_old-p*fn; mom_new.setEnergy( pf4 );
  setMomentum( mom_new );
}

void
Particle::Momentum::rotatePhi( double phi, double sign )
{
  const double px = px_*cos( phi )+py_*sin( phi )*sign,
               py =-px_*sin( phi )+py_*cos( phi )*sign;
  px_ = px;
  py_ = py;
}

void
Particle::Momentum::rotateThetaPhi( double theta, double phi )
{
  double rotmtx[3][3], mom[3]; //FIXME check this! cos(phi)->-sin(phi) & sin(phi)->cos(phi) --> phi->phi+pi/2 ?
  rotmtx[0][0] = -sin( phi ); rotmtx[0][1] = -cos( theta )*cos( phi ); rotmtx[0][2] =  sin( theta )*cos( phi );
  rotmtx[1][0] =  cos( phi ); rotmtx[1][1] = -cos( theta )*sin( phi ); rotmtx[1][2] =  sin( theta )*sin( phi );
  rotmtx[2][0] =  0.;         rotmtx[2][1] =  sin( theta );            rotmtx[2][2] =  cos( theta );

  for (int i=0; i<3; i++) {
    mom[i] = 0.;
    for (int j=0; j<3; j++) {
      mom[i] += rotmtx[i][j]*p( j );
    }
  }

  setP( mom[0], mom[1], mom[2] );
}

std::ostream&
operator<<( std::ostream& os, const Particle::Momentum& mom )
{
  os << "(E, p) = (" << mom.energy_ << ", " << mom.px_ << ", " << mom.py_ << ", " << mom.pz_ << ")";
  return os;
}

