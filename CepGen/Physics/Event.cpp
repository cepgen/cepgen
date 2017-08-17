#include "Event.h"

namespace CepGen
{
  Event::Event() :
    num_hadronisation_trials( 0 ),
    time_generation( -1. ), time_total( -1. )
  {}

  Event::~Event()
  {}

  Event&
  Event::operator=( const Event &ev_ )
  {
    particles_ = ev_.particles_;
    time_generation = ev_.time_generation;
    time_total = ev_.time_total;
    num_hadronisation_trials = ev_.num_hadronisation_trials;
    return *this;
  }

  void
  Event::clear()
  {
    particles_.clear();
    time_generation = -1.;
    time_total = -1.;
  }

  void
  Event::init()
  {
    last_particle_ = particles_.end();
  }

  void
  Event::restore()
  {
    //--- remove all particles after the primordial event block
    particles_.erase( last_particle_, particles_.end() );
  }

  Particles&
  Event::getByRole( const Particle::Role& role )
  {
    //--- retrieve all particles with a given role
    return particles_[role];
  }

  Particle&
  Event::getOneByRole( const Particle::Role& role )
  {
    //--- retrieve the first particle a the given role
    Particles& parts_by_role = getByRole( role );
    if ( parts_by_role.size() > 1 ) {
      InWarning( Form( "More than one particle with role %d: %d particles", (int)role, parts_by_role.size() ) );
    }
    return *parts_by_role.begin();
  }

  Particle&
  Event::getById( int id )
  {
    for ( ParticlesMap::iterator out=particles_.begin(); out!=particles_.end(); out++ ) {
      for ( Particles::iterator part=out->second.begin(); part!=out->second.end(); part++ ) {
        if ( part->id == id ) return *part;
      }
    }
    throw Exception( __PRETTY_FUNCTION__, Form( "Failed to retrieve the particle with id=%d", id ), FatalError );
  }

  const Particle&
  Event::getConstById( int id ) const
  {
    for ( ParticlesMap::const_iterator out=particles_.begin(); out!=particles_.end(); out++ ) {
      for ( Particles::const_iterator part=out->second.begin(); part!=out->second.end(); part++ ) {
        if ( part->id == id ) return *part;
      }
    }
    throw Exception( __PRETTY_FUNCTION__, Form( "Failed to retrieve the particle with id=%d", id ), FatalError );
  }

  Particles
  Event::getByIds( const ParticlesIds& ids ) const
  {
    Particles out;
    for ( ParticlesIds::const_iterator id=ids.begin(); id!=ids.end(); id++ ) {
      out.emplace_back( getConstById( *id ) );
    }
    return out;
  }

  Particles
  Event::mothers( const Particle& part )
  {
    return getByIds( part.mothersIds() );
  }

  Particles
  Event::daughters( const Particle& part )
  {
    return getByIds( part.daughters() );
  }

  ParticleRoles
  Event::roles() const
  {
    ParticleRoles out;
    ParticlesMap::const_iterator it, end;
    for ( it=particles_.begin(), end=particles_.end(); it!=end; it=particles_.upper_bound( it->first ) ) {
      out.emplace_back( it->first );
    }
    return out;
  }

  int
  Event::addParticle( Particle part, bool replace )
  {
    DebuggingInsideLoop( Form( "Particle with PDGid = %d has role %d", part.pdgId(), part.role ) );
    if ( part.role <= 0 ) return -1;

    Particles& part_with_same_role = getByRole( part.role );
    part.id = particles().size(); //FIXME is there any better way of introducing this id ?
    if ( !replace ) part_with_same_role.emplace_back( part );
    else part_with_same_role = Particles( 1, part );
    return 1;
  }

  int
  Event::addParticle( const Particle::Role& role, bool replace )
  {
    if ( role <= 0 ) return -1;
    return addParticle( Particle( role ), replace );
  }

  Particles
  Event::particles() const
  {
    Particles out;
    for ( ParticlesMap::const_iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
      for ( Particles::const_iterator part=it->second.begin(); part!=it->second.end(); part++ ) {
        out.emplace_back( *part );
      }
    }
    std::sort( out.begin(), out.end(), compareParticle );
    return out;
  }

  Particles
  Event::stableParticles() const
  {
    Particles out;
    for ( ParticlesMap::const_iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
      for ( Particles::const_iterator part=it->second.begin(); part!=it->second.end(); part++ ) {
        if ( part->status==Particle::Undefined
          || part->status==Particle::FinalState ) {
          out.emplace_back( *part );
        }
      }
    }
    std::sort( out.begin(), out.end(), compareParticle );
    return out;
  }

  void
  Event::dump( bool stable ) const
  {
    double pxtot, pytot, pztot, etot;
    std::ostringstream os;

    pxtot = pytot = pztot = etot = 0.;
    const Particles parts = particles();
    for ( Particles::const_iterator part_ref=parts.begin(); part_ref!=parts.end(); part_ref++ ) {
      const Particle& p = *part_ref;
      if ( stable && p.status != Particle::FinalState ) continue;

      os << Form( "\n %2d\t%+6d", p.id, p.integerPdgId() );
      if ( p.name!="" ) os << Form( "%6s", p.name.c_str() );
      else              os << "\t";
      os << "\t";
      if ( p.charge!=999. ) os << Form( "%6.2f\t", p.charge );
      else                  os << "\t";
      os << Form( "%4d\t%6d\t", p.role, p.status );
      if ( !p.mothersIds().empty() )
        os << Form( "%2d(%2d)", *( p.mothersIds().begin() ), getConstById( *( p.mothersIds().begin() ) ).role );
      else
        os << "      ";
      os << Form( "% 9.3e % 9.3e % 9.3e % 9.3e", p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.energy() );
      if ( p.status == Particle::Undefined
        || p.status == Particle::FinalState
        || p.status == Particle::Undecayed ) {
        const int sign = ( p.status == Particle::Undefined ) ? -1 : 1;
        pxtot += sign*p.momentum().px();
        pytot += sign*p.momentum().py();
        pztot += sign*p.momentum().pz();
        etot += sign*p.energy();
      }
    }
    //--- set a threshold to the computation precision
    if ( fabs( pxtot ) < 1.e-12 ) pxtot = 0.;
    if ( fabs( pytot ) < 1.e-12 ) pytot = 0.;
    if ( fabs( pztot ) < 1.e-12 ) pztot = 0.;
    if ( fabs(  etot ) < 1.e-12 ) etot = 0.;
    //
    Information( Form( "Dump of event content:\n"
    "Part.\tPDG id\t\tCharge\tRole\tStatus\tMother\t\t4-Momentum (GeV)\n"
    "----\t------\t\t------\t----\t------\t------\t-------------------------------------"
    "%s\n"
    "---------------------------------------------------------------------------------------------\n"
    "Total:\t\t\t\t\t\t      % 9.4f % 9.4f % 9.4f % 9.4f", os.str().c_str(), pxtot, pytot, pztot, etot ) );
  }
}
