#include "Particle.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Physics/Constants.h"

namespace CepGen
{
  Particle::Particle() :
    id_( -1 ), charge_sign_( 1 ),
    mass_( -1. ), helicity_( 0. ),
    role_( UnknownRole ), status_( Undefined ), pdg_id_( invalidParticle )
  {}

  Particle::Particle( Role role, ParticleCode pdgId, Status st ) :
    id_( -1 ), charge_sign_( 1 ),
    mass_( -1. ), helicity_( 0. ),
    role_( role ), status_( st ), pdg_id_( pdgId )
  {
    if ( pdg_id_!=invalidParticle ) {
      computeMass();
    }
  }

  Particle::Particle( const Particle& part ) :
    id_( part.id_ ), charge_sign_( part.charge_sign_ ),
    momentum_( part.momentum_ ), mass_( part.mass_ ), helicity_( part.helicity_ ),
    role_( part.role_ ), status_( part.status_ ),
    mothers_( part.mothers_ ), daughters_( part.daughters_ ),
    pdg_id_( part.pdg_id_ )
  {}

  bool
  Particle::operator<( const Particle& rhs ) const
  {
    return ( id_ >= 0 && rhs.id_ > 0 && id_ < rhs.id_ );
  }

  double
  Particle::thetaToEta( double theta )
  {
    return -log( tan( 0.5 * theta * M_PI/180. ) );
  }

  double
  Particle::etaToTheta( double eta )
  {
    return 2.*atan( exp( -eta ) )*180. / M_PI;
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
      mass_ = ParticleProperties::mass( pdg_id_ );
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

    DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is the new mother of %2d (pdgId=%4d)",
                               part.id(), part.pdgId(), id_, pdg_id_ ) );

    part.addDaughter( *this );
  }

  void
  Particle::addDaughter( Particle& part )
  {
    std::pair<ParticlesIds::iterator,bool> ret = daughters_.insert( part.id() );

    if ( Logger::get().level >= Logger::DebugInsideLoop ) {
      std::ostringstream os;
      for ( ParticlesIds::const_iterator it = daughters_.begin(); it != daughters_.end(); ++it ) {
        os << Form( "\n\t * id=%d", *it );
      }
      DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) has now %2d daughter(s):"
                                 "%s", role_, pdg_id_, numDaughters(), os.str().c_str() ) );
    }

    if ( ret.second ) {
      DebuggingInsideLoop( Form( "Particle %2d (pdgId=%4d) is a new daughter of %2d (pdgId=%4d)",
                                 part.role(), part.pdgId(), role_, pdg_id_ ) );

      if ( part.mothers().find( id_ ) == part.mothers().end() ) part.addMother( *this );
    }
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

  double
  Particle::energy() const
  {
    return ( momentum_.energy() < 0. ) ? sqrt( mass2()+momentum_.p2() ) : momentum_.energy();
  }

  void
  Particle::setEnergy( double e )
  {
    if ( e < 0. && mass_ >= 0. ) e = sqrt( mass2()+momentum_.p2() );
    momentum_.setEnergy( e );
  }

  void
  Particle::setPdgId( const ParticleCode& pdg, short ch )
  {
    pdg_id_ = pdg;
    charge_sign_ = ch;
  }

  int
  Particle::integerPdgId() const
  {
    const float ch = ParticleProperties::charge( pdg_id_ );
    if ( ch == 0 ) return static_cast<int>( pdg_id_ );
    return static_cast<int>( pdg_id_ ) * charge_sign_ * ( ch/fabs( ch ) );
  }

  void
  Particle::dump() const
  {
    std::ostringstream osm, osd, os;
    if ( !primary() ) {
      osm << ": mother(s): ";
      for ( ParticlesIds::const_iterator m = mothers_.begin(); m != mothers_.end(); ++m ) {
        if ( m!=mothers_.begin() ) osm << ", ";
        osm << ( *m );
      }
    }
    const ParticlesIds daugh = daughters();
    if ( daugh.size() != 0 ) {
      osd << ": id = ";
      for ( ParticlesIds::const_iterator it = daugh.begin(); it != daugh.end(); ++it ) {
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
    if ( mom.energy() == m )
      return *this;

    double pf4 = 0.;
    for ( unsigned int i = 0; i < 4; ++i )
      pf4 += momentum_[i]*mom[i];
    pf4 /= m;
    const double fn = ( pf4+energy() )/( momentum_.energy()+m );
    for ( unsigned int i = 0; i < 3; ++i )
      momentum_.setP( i, momentum_[i]+fn*mom[i] );

    return *this;
  }

  std::vector<double>
  Particle::lorentzBoost( const Particle::Momentum& mom )
  {
    std::vector<double> out( 3, 0. );

    const double p2 = mom.p2();
    const double gamma = 1./sqrt( 1.-p2 );
    double bp = 0.;
    for ( unsigned int i = 0; i < 3; ++i )
      bp += mom[i]*momentum_[i];

    const double gamma2 = ( p2 > 0. ) ? ( gamma-1. )/p2 : 0.;

    for ( unsigned int i = 0; i < 3; ++i ) {
      out[i] = momentum_[i] + gamma2*bp*mom[i]+gamma*mom[i]*energy();
    }
    return out;
  }

  double
  Particle::etaToY( double eta_, double m_, double pt_ )
  {
    const double mt = m_*m_ + pt_*pt_;
    return asinh( sqrt( ( ( ( mt*mt-m_*m_ )*cosh( 2.*eta_ ) + m_*m_ )/ mt*mt - 1. ) / 2. ) );
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
      case Particle::Intermediate:  return os << "hard.pr.";
      case Particle::CentralSystem: return os << "central";
    }
    return os;
  }

  double
  CMEnergy( const Particle& p1, const Particle& p2 )
  {
    if ( p1.mass()*p2.mass() < 0. ) return 0.;
    if ( p1.energy()*p2.energy() < 0. ) return 0.;
    return sqrt( p1.mass2()+p2.mass2() + 2.*p1.energy()*p2.energy() - 2.*( p1.momentum()*p2.momentum() ) );
  }

  double
  CMEnergy( const Particle::Momentum& m1, const Particle::Momentum& m2 )
  {
    if ( m1.mass()*m2.mass() < 0. ) return 0.;
    if ( m1.energy()*m2.energy() < 0. ) return 0.;
    return sqrt( m1.mass2()+m2.mass2() + 2.*m1.energy()*m2.energy() - 2.*( m1*m2 ) );
  }
}
