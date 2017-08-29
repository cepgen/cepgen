#include "Particle.h"

namespace CepGen
{
  Particle::Particle() :
    id_( -1 ), charge_( 1. ),
    mass_( -1. ), helicity_( 0. ),
    role_( UnknownRole ), status_( Undefined ), pdg_id_( invalidParticle ), is_primary_( true )
  {}

  Particle::Particle( Role role, ParticleCode pdgId ) :
    id_( -1 ), charge_( 1. ),
    mass_( -1. ), helicity_( 0. ),
    role_( role ), status_( Undefined ), pdg_id_( pdgId ), is_primary_( true )
  {
    if ( pdg_id_!=invalidParticle ) {
      computeMass();
    }
  }

  Particle::Particle( const Particle& part ) :
    id_( part.id_ ), charge_( part.charge_ ),
    momentum_( part.momentum_ ), mass_( part.mass_ ), helicity_( part.helicity_ ),
    role_( part.role_ ), status_( part.status_ ),
    mothers_( part.mothers_ ), daughters_( part.daughters_ ),
    pdg_id_( part.pdg_id_ ), is_primary_( part.is_primary_ )
  {}

  bool
  Particle::operator<( const Particle& rhs ) const
  {
    return ( id_ < rhs.id_ );
  }

  bool
  Particle::valid()
  {
    if ( pdg_id_ == invalidParticle ) return false;
    if ( momentum_.p() == 0. && mass() == 0. ) return false;
    return true;
  }

  void
  Particle::computeMass( bool off_shell )
  {
    if ( !off_shell && pdg_id_ != invalidParticle ) { // retrieve the mass from the on-shell particle's properties
      mass_ = massFromPDGId( pdg_id_ );
    }
    else if ( momentum_.energy() >= 0. ) {
      mass_ = sqrt( energy2() - momentum_.p2() );
    }

    //--- finish by setting the energy accordingly
    if ( momentum_.energy() < 0. ) { // invalid energy
      momentum_.setEnergy( sqrt( momentum_.p2() + mass2() ) );
    }
  }

  void
  Particle::setMass( double m )
  {
    if ( m >= 0. ) mass_ = m;
    else computeMass();
  }

  void
  Particle::addMother( Particle& part )
  {
    mothers_.insert( part.id() );
    is_primary_ = false;

    DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is the new mother of %2d (pdgId=%4d)",
                               part.id()+1, part.pdgId(), id_+1, pdg_id_ ) );

    part.addDaughter( *this );
  }

  bool
  Particle::addDaughter( Particle& part )
  {
    std::pair<ParticlesIds::iterator,bool> ret = daughters_.insert( part.id() );

    if ( Logger::get().level >= Logger::DebugInsideLoop ) {
      std::ostringstream os;
      for ( ParticlesIds::const_iterator it=daughters_.begin(); it!=daughters_.end(); it++) {
        os << Form("\n\t * id=%d", *it);
      }
      DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) has now %2d daughter(s):"
                                 "%s", role_, pdg_id_, numDaughters(), os.str().c_str() ) );
    }

    if ( ret.second ) {
      DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is a new daughter of %2d (pdgId=%4d)",
                                 part.role(), part.pdgId(), role_, pdg_id_ ) );

      if ( !part.primary() && part.mothersIds().empty() ) {
        part.addMother( *this );
      }
    }

    return ret.second;
  }

  void
  Particle::setMomentum( const Momentum& mom, bool offshell )
  {
    momentum_ = mom;
    if ( !offshell && mom.mass() > 0. ) mass_ = momentum_.mass();
    else computeMass();
  }

  void
  Particle::setMomentum( double px, double py, double pz )
  {
    momentum_.setP( px, py, pz );
    setEnergy();
  }

  void
  Particle::setMomentum( double px, double py, double pz, double e )
  {
    setMomentum( px, py, pz );
    if ( fabs( e-momentum_.energy() )>1.e-6 ) { // more than 1 eV difference
      InError( Form( "Energy difference: %.5e", e-momentum_.energy() ) );
      return;
    }
  }

  void
  Particle::setEnergy( double e )
  {
    if ( e < 0. && mass_ >= 0. ) e = sqrt( mass2()+momentum_.p2() );
    momentum_.setEnergy( e );
  }

  void
  Particle::setPdgId( const ParticleCode& pdg, float ch )
  {
    pdg_id_ = pdg;
    if ( ch == -999. ) charge_ = 0.;
    else charge_ = ch;
  }

  int
  Particle::integerPdgId() const
  {
    const int pdg = static_cast<int>( pdg_id_ );
    //--- leptons
    if ( charge_ != 0 && pdg > 10 && pdg < 16 && pdg%2 != 0 ) return -charge_*pdg;
    return pdg;
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
    if ( os.str() == " ()" ) os.str("");
    Information( Form(
      "Dumping a particle with id=%3d, role=%3d, status=% 3d\n\t"
      "PDG Id:%4d%s, mass = %5.4f GeV\n\t"
      "(E,P) = (%4.2f, %4.2f, %4.2f, %4.2f) GeV\t"
      "(|P| = p = %4.2f GeV)\n\t"
      " Pt = %5.4f GeV, eta = %4.3f, phi = % 4.3f\n\t"
      "Primary? %s%s\n\t"
      "%d daughter(s)%s",
      id_, role_, status_, pdg_id_, os.str().c_str(),
      mass(), energy(), momentum_.px(), momentum_.py(), momentum_.pz(),
      momentum_.p(), momentum_.pt(), momentum_.eta(), momentum_.phi(),
      yesno( primary() ), osm.str().c_str(), numDaughters(), osd.str().c_str() )
    );
  }

  Particle&
  Particle::lorentzBoost( double m, const Particle::Momentum& mom )
  {
    double pf4, fn;

    if ( mom.energy() != m ) {
      pf4 = 0.;
      for ( unsigned int i=0; i<4; i++ ) {
        pf4 += momentum_[i]*mom[i];
      }
      pf4 /= m;
      fn = ( pf4+energy() )/( momentum_.energy()+m );
      for ( unsigned int i=0; i<3; i++ ) {
        momentum_.setP( i, momentum_[i]+fn*mom[i] );
      }
    }
    return *this;
  }

  std::vector<double>
  Particle::lorentzBoost( const Particle::Momentum& mom )
  {
    double p2, gamma, bp, gamma2;
    std::vector<double> out( 3, 0. );

    p2 = mom.p2();
    gamma = 1./sqrt( 1.-p2 );
    bp = 0.;
    for ( unsigned int i=0; i<3; i++ ) bp+= mom[i]*momentum_[i];

    if ( p2>0. ) gamma2 = (gamma-1.)/p2;
    else gamma2 = 0.;

    for ( unsigned int i=0; i<3; i++ ) {
      out[i] = momentum_[i] + gamma2*bp*mom[i]+gamma*mom[i]*energy();
    }
    return out;
  }

  double
  Particle::massFromPDGId( const Particle::ParticleCode& pdgId_ )
  {
    switch ( pdgId_ ) {
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
      case Particle::dQuark:          return os << "d quark";
      case Particle::uQuark:          return os << "u quark";
      case Particle::Electron:        return os << "electron";
      case Particle::ElectronNeutrino: return os << "electron neutrino";
      case Particle::Muon:            return os << "muon";
      case Particle::MuonNeutrino:    return os << "muon neutrino";
      case Particle::Tau:             return os << "tau";
      case Particle::TauNeutrino:     return os << "tau neutrino";
      case Particle::Gluon:           return os << "gluon";
      case Particle::Photon:          return os << "photon";
      case Particle::Z:               return os << "Z";
      case Particle::WPlus:           return os << "W+";
      case Particle::PiPlus:          return os << "pi+";
      case Particle::PiZero:          return os << "pi0";
      case Particle::Rho770_0:        return os << "rho(770)0";
      case Particle::Rho1450_0:       return os << "rho(1450)0";
      case Particle::Rho1700_0:       return os << "rho(1700)0";
      case Particle::h1380_1:         return os << "h(1380)1";
      case Particle::Omega782:        return os << "omega(782)";
      case Particle::JPsi:            return os << "J/Psi";
      case Particle::Phi1680:         return os << "phi(1680)";
      case Particle::Upsilon1S:       return os << "Upsilon(1S)";
      case Particle::Upsilon2S:       return os << "Upsilon(2S)";
      case Particle::Upsilon3S:       return os << "Upsilon(3S)";;
      case Particle::ud0Diquark:      return os << "(ud)0 di-quark";
      case Particle::ud1Diquark:      return os << "(ud)1 di-quark";
      case Particle::uu1Diquark:      return os << "(uu)1 di-quark";
      case Particle::Proton:          return os << "proton";
      case Particle::Neutron:         return os << "neutron";
      case Particle::Pomeron:         return os << "pomeron";
      case Particle::Reggeon:         return os << "reggeon";
      case Particle::invalidParticle: return os << "[...]";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Particle::Role& rl )
  {
    switch ( rl ) {
      case Particle::UnknownRole:   return os << "unknown";
      case Particle::IncomingBeam1: return os << "in.b.1";
      case Particle::IncomingBeam2: return os << "in.b.2";
      case Particle::OutgoingBeam1: return os << "out.b.1";
      case Particle::OutgoingBeam2: return os << "out.b.2";
      case Particle::Parton1:       return os << "parton1";
      case Particle::Parton2:       return os << "parton2";
      case Particle::Parton3:       return os << "parton3";
      case Particle::Intermediate:  return os << "partons";
      case Particle::CentralSystem: return os << "central";
    }
    return os;
  }
}
