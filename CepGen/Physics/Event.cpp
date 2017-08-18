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
      FatalError( Form( "More than one particle with role %d: %d particles", (int)role, parts_by_role.size() ) );
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

  void
  Event::addParticle( Particle part, bool replace )
  {
    DebuggingInsideLoop( Form( "Particle with PDGid = %d has role %d", part.pdgId(), part.role ) );
    if ( part.role <= 0 ) FatalError( Form( "Trying to add a particle with role=%d", (int)part.role ) );

    //--- retrieve the list of particles with the same role
    Particles& part_with_same_role = getByRole( part.role );

    //--- specify the id
    if ( part_with_same_role.empty() && part.id < 0 ) part.id = particles().size(); //FIXME is there any better way of introducing this id ?
    if ( replace && !part_with_same_role.empty() ) part.id = part_with_same_role[0].id;

    //--- add the particle to the collection
    if ( replace ) part_with_same_role = Particles( 1, part ); // generate a vector containing only this particle
    else part_with_same_role.emplace_back( part );

  }

  void
  Event::addParticle( const Particle::Role& role, bool replace )
  {
    addParticle( Particle( role ), replace );
  }

  const Particles
  Event::particles() const
  {
    Particles out;
    for ( ParticlesMap::const_iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
      out.insert( out.end(), it->second.begin(), it->second.end() );
    }
    std::cout<<"before:"<<std::endl;for (const auto& p:out) std::cout<<p.id<<" "; std::cout<<std::endl;
    std::sort( out.begin(), out.end() );
    std::cout<<"after:"<<std::endl;for (const auto& p:out) std::cout<<p.id <<" "; std::cout<<std::endl;
    std::cout<<"sorted? "<< std::boolalpha << std::is_sorted(out.begin(),out.end())<<std::endl;
    return out;
  }

  const Particles
  Event::stableParticles() const
  {
    Particles out;
    for ( ParticlesMap::const_iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
      for ( Particles::const_iterator part=it->second.begin(); part!=it->second.end(); part++ ) {
        if ( part->status == Particle::Undefined || part->status == Particle::FinalState ) {
          out.emplace_back( *part );
        }
      }
    }
    std::sort( out.begin(), out.end() );
    return out;
  }

  void
  Event::dump( std::ostream& out, bool stable ) const
  {
    Particles parts = particles();

    std::ostringstream os;

    double pxtot = 0., pytot = 0., pztot = 0., etot = 0.;
    for ( Particles::const_iterator part_ref=parts.begin(); part_ref!=parts.end(); part_ref++ ) {
      const Particle& p = *part_ref;
      if ( stable && p.status != Particle::FinalState ) continue;

      std::ostringstream oss; oss << p.pdgId();
      os << Form( "\n %2d\t%+6d%8s", p.id, p.integerPdgId(), oss.str().c_str() );
      os << "\t";
      if ( p.charge!=999. ) os << Form( "%6.2f\t", p.charge );
      else                  os << "\t";
      os << Form( "%4d\t%6d\t", p.role, p.status );
      if ( !p.mothersIds().empty() )
        os << Form( "%2d (%2d)", *( p.mothersIds().begin() ), getConstById( *( p.mothersIds().begin() ) ).role );
      else
        os << "       ";
      os << Form( "% 9.6e % 9.6e % 9.6e % 9.6e % 9.5e", p.momentum().px(), p.momentum().py(), p.momentum().pz(), p.energy(), p.mass() );
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
    "Part.\tPDG id\t\tCharge\tRole\tStatus\tMother\t\t\t\t4-Momentum (GeV)\t\tMass (GeV)\n"
    "----\t------\t\t------\t----\t------\t------\t------------------------------------------------------  -----------"
    "%s\n"
    "---------------------------------------------------------------------------------------------------------------------------\n"
    "\t\t\t\t\t\tTotal: % 9.6e % 9.6e % 9.6e % 9.6e", os.str().c_str(), pxtot, pytot, pztot, etot ) );
  }
}
