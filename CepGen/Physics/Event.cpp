#include "Event.h"

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
  particles_.erase( last_particle_, particles_.end() );
}

ParticlesRef
Event::getByRole( const Particle::Role& role_ )
{
  ParticlesRef out;
  std::pair<ParticlesMap::iterator,ParticlesMap::iterator> ret = particles_.equal_range( role_ );

  unsigned int i = 0;
  for ( ParticlesMap::iterator it=ret.first; it!=ret.second && i<100; it++ ) {
    out.push_back( &( it->second ) );
    i++;
  }
  return out;
}

Particle*
Event::getById( int id_ )
{
  for ( ParticlesMap::iterator out=particles_.begin(); out!=particles_.end(); out++ ) {
    if ( out->second.id==id_ ) return &out->second;
  }
  return 0;
}

const Particle
Event::getConstById( int id_ ) const
{
  for ( ParticlesMap::const_iterator out=particles_.begin(); out!=particles_.end(); out++ ) {
    if ( out->second.id==id_ ) return static_cast<Particle>( out->second );
  }
  return Particle();
}

ParticleRoles
Event::roles() const
{
  ParticleRoles out;
  ParticlesMap::const_iterator it, end;
  for ( it=particles_.begin(), end=particles_.end(); it!=end; it=particles_.upper_bound( it->first ) ) {
    out.push_back( it->first );
  }
  return out;
}

int
Event::addParticle( Particle part_, bool replace_ )
{
  DebuggingInsideLoop( Form( "Particle with PDGid = %d has role %d", part_.pdgId(), part_.role ) );
  if ( part_.role<=0 ) return -1;

  ParticlesRef part_with_same_role = getByRole( part_.role );
  part_.id = particles_.size(); //FIXME is there any better way of introducing this id ?
  if ( replace_ and part_with_same_role.size()!=0 ) {
    part_with_same_role.at( 0 ) = &part_;
    return 0;
  }
  particles_.insert( std::pair<Particle::Role,Particle>( part_.role, part_ ) );
  return 1;
}

int
Event::addParticle( const Particle::Role& role_, bool replace_ )
{
  if ( role_<=0 ) return -1;
  return addParticle( Particle( role_ ), replace_ );
}

ParticlesRef
Event::particles()
{
  ParticlesRef out;
  for ( ParticlesMap::iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
    out.push_back( &it->second );
  }
  std::sort( out.begin(), out.end(), compareParticlePtrs );
  return out;
}

Particles
Event::constParticles() const
{
  Particles out;
  for ( ParticlesMap::const_iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
    out.push_back( static_cast<Particle>( it->second ) );
  }
  std::sort( out.begin(), out.end(), compareParticle );
  return out;
}

ConstParticlesRef
Event::constParticlesRef() const
{
  ConstParticlesRef out;
  for ( ParticlesMap::const_iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
    out.push_back( &it->second );
  }
  std::sort( out.begin(), out.end(), compareParticlePtrs );
  return out;
}

ParticlesRef
Event::stableParticles()
{
  ParticlesRef out;
  for ( ParticlesMap::iterator it=particles_.begin(); it!=particles_.end(); it++ ) {
    if ( it->second.status==Particle::Undefined
      or it->second.status==Particle::FinalState ) {
      out.push_back( &it->second );
    }
  }
  std::sort( out.begin(), out.end() );
  return out;  
}

void
Event::dump( bool stable ) const
{
  double pxtot, pytot, pztot, etot;
  std::ostringstream os;

  pxtot = pytot = pztot = etot = 0.;
  const ConstParticlesRef particles = constParticlesRef();
  for ( ConstParticlesRef::const_iterator part_ref=particles.begin(); part_ref!=particles.end(); part_ref++ ) {
    const Particle* p = ( *part_ref );
    if ( stable and p->status!=Particle::FinalState ) continue;

    os << Form( "\n %2d\t%+6d", p->id, p->integerPdgId() );
    if ( p->name!="" ) os << Form( "%6s", p->name.c_str() );
    //else               os << std::setw(6) << Particle::ParticleCode( abs( p->pdgId() ) );
    else               os << "\t";
    os << "\t";
    if ( p->charge!=999. ) os << Form( "%6.2f\t", p->charge );
    else                   os << "\t";
    os << Form( "%4d\t%6d\t", p->role, p->status );
    if ( p->mothersIds().size()>0 ) 
      os << Form( "%2d(%2d)", *( p->mothersIds().begin() ), getConstById( *( p->mothersIds().begin() ) ).role );
    else
      os << "      ";
    os << Form( "% 9.3f % 9.3f % 9.3f % 9.3f", p->momentum().px(), p->momentum().py(), p->momentum().pz(), p->energy() );
    if ( p->status==Particle::Undefined
      or p->status==Particle::FinalState
      or p->status==Particle::Undecayed ) {
      const int sign = ( p->status==Particle::Undefined ) ? -1 : 1;
      pxtot += sign*p->momentum().px();
      pytot += sign*p->momentum().py();
      pztot += sign*p->momentum().pz();
      etot += sign*p->energy();
    }
  }
  // We set a threshold to the computation precision
  if ( fabs( pxtot )<1.e-12 ) pxtot = 0.;
  if ( fabs( pytot )<1.e-12 ) pytot = 0.;
  if ( fabs( pztot )<1.e-12 ) pztot = 0.;
  if ( fabs(  etot )<1.e-12 ) etot = 0.;
  //
  Information( Form( "Dump of event content:\n"
  "Part.\tPDG id\t\tCharge\tRole\tStatus\tMother\t\t4-Momentum (GeV)\n"
  "----\t------\t\t------\t----\t------\t------\t-------------------------------------"
  "%s\n"
  "---------------------------------------------------------------------------------------------\n"
  "Total:\t\t\t\t\t\t      % 9.4f % 9.4f % 9.4f % 9.4f", os.str().c_str(), pxtot, pytot, pztot, etot ) );
}
