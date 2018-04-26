#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Physics/Constants.h"

namespace CepGen
{
  Particle::Particle() :
    id_( -1 ), charge_sign_( 1 ),
    mass_( -1. ), helicity_( 0. ),
    role_( UnknownRole ), status_( Undefined ), pdg_id_( PDG::invalid )
  {}

  Particle::Particle( Role role, PDG pdgId, Status st ) :
    id_( -1 ), charge_sign_( 1 ),
    mass_( -1. ), helicity_( 0. ),
    role_( role ), status_( st ), pdg_id_( pdgId )
  {
    if ( pdg_id_ != PDG::invalid ) {
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
    if ( pdg_id_ == PDG::invalid ) return false;
    if ( momentum_.p() == 0. && mass() == 0. ) return false;
    return true;
  }

  void
  Particle::computeMass( bool off_shell )
  {
    if ( !off_shell && pdg_id_ != PDG::invalid ) { // retrieve the mass from the on-shell particle's properties
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

    CG_DEBUG_LOOP( "Particle" )
      <<  "Particle " << id() << " (pdgId=" << part.integerPdgId() << ") "
      << "is the new mother of " << id_ << " (pdgId=" << (int)pdg_id_ << ").";

    part.addDaughter( *this );
  }

  void
  Particle::addDaughter( Particle& part )
  {
    const auto ret = daughters_.insert( part.id() );

    if ( Logger::get().level >= Logger::Level::DebugInsideLoop ) {
      std::ostringstream os;
      for ( const auto& daugh : daughters_ )
        os << Form( "\n\t * id=%d", daugh );
      CG_DEBUG_LOOP( "Particle" )
        << "Particle " << role_ << " (pdgId=" << (int)pdg_id_ << ") "
        << "has now " << daughters_.size() << " daughter(s):"
        << os.str();
    }

    if ( ret.second ) {
      CG_DEBUG_LOOP( "Particle" )
        << "Particle " << part.role() << " (pdgId=" << part.integerPdgId() << ") "
        << "is a new daughter of " << role_ << " (pdgId=" << (int)pdg_id_ << "%4d).";

      if ( part.mothers().find( id_ ) == part.mothers().end() )
        part.addMother( *this );
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
      CG_ERROR( Form( "Energy difference: %.5e", e-momentum_.energy() ) );
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
  Particle::setPdgId( const PDG& pdg, short ch )
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
    std::ostringstream osm, osd;
    if ( !primary() ) {
      osm << ": mother(s): ";
      unsigned short i = 0;
      for ( const auto& moth : mothers_ ) {
        osm << ( i > 0 ? ", " : "" ) << moth;
        ++i;
      }
    }
    const ParticlesIds daughters_list = daughters();
    if ( daughters_list.size() > 0 ) {
      osd << ": id = ";
      unsigned short i = 0;
      for ( const auto& daugh : daughters_list ) {
        osm << ( i > 0 ? ", " : "" ) << daugh;
        ++i;
      }
    }
    CG_INFO( "Particle" )
      << "Dumping a particle with id=" << id_ << ", role=" << role_ << ", status=" << status_ << "\n\t"
      << "Particle id: " << integerPdgId() << " (" << pdg_id_ << "), mass = " << mass() << " GeV\n\t"
      << "Momentum: " << momentum_ << " GeV\t" << "(|P| = p = " << momentum_.p() << " GeV)\n\t"
      << " p⟂ = " << momentum_.pt() << " GeV, eta = " << momentum_.eta() << ", phi = " << momentum_.phi() << "\n\t"
      << "Primary? " << yesno( primary() ) << osm.str() << "\n\t"
      << numDaughters() << " daughter(s)" << osd.str();
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
