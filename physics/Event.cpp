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
  fParticles = ev_.fParticles;
  time_generation = ev_.time_generation;
  time_total = ev_.time_total;
  num_hadronisation_trials = ev_.num_hadronisation_trials;
  return *this;
}

void
Event::clear()
{
  fParticles.clear();
  time_generation = -1.;
  time_total = -1.;
}

void
Event::Init()
{
  fLastParticle = fParticles.end();
}

void
Event::Restore()
{
  fParticles.erase( fLastParticle, fParticles.end() );
}

ParticlesRef
Event::GetByRole( const Particle::Role& role_ )
{
  ParticlesRef out;
  std::pair<ParticlesMap::iterator,ParticlesMap::iterator> ret = fParticles.equal_range( role_ );

  unsigned int i = 0;
  for ( ParticlesMap::iterator it=ret.first; it!=ret.second && i<100; it++ ) {
    out.push_back( &( it->second ) );
    i++;
  }
  return out;
}

Particle*
Event::GetById( int id_ )
{
  for ( ParticlesMap::iterator out=fParticles.begin(); out!=fParticles.end(); out++ ) {
    if ( out->second.id==id_ ) return &out->second;
  }
  return 0;
}

const Particle
Event::GetConstById( int id_ ) const
{
  for ( ParticlesMap::const_iterator out=fParticles.begin(); out!=fParticles.end(); out++ ) {
    if ( out->second.id==id_ ) return static_cast<Particle>( out->second );
  }
  return Particle();
}

ParticleRoles
Event::GetRoles() const
{
  ParticleRoles out;
  ParticlesMap::const_iterator it, end;
  for ( it=fParticles.begin(), end=fParticles.end(); it!=end; it=fParticles.upper_bound( it->first ) ) {
    out.push_back( it->first );
  }
  return out;
}

int
Event::AddParticle( Particle part_, bool replace_ )
{
  DebuggingInsideLoop( Form( "Particle with PDGid = %d has role %d", part_.GetPDGId(), part_.role ) );
  if ( part_.role<=0 ) return -1;

  ParticlesRef part_with_same_role = GetByRole( part_.role );
  part_.id = fParticles.size(); //FIXME is there any better way of introducing this id ?
  if ( replace_ and part_with_same_role.size()!=0 ) {
    part_with_same_role.at( 0 ) = &part_;
    return 0;
  }
  fParticles.insert( std::pair<Particle::Role,Particle>( part_.role, part_ ) );
  return 1;
}

int
Event::AddParticle( const Particle::Role& role_, bool replace_ )
{
  if ( role_<=0 ) return -1;

  np = new Particle();
  np->role = role_;
  int out = AddParticle( *np, replace_ );

  delete np;
  return out;
}

ParticlesRef
Event::GetParticles()
{
  ParticlesRef out;
  for ( ParticlesMap::iterator it=fParticles.begin(); it!=fParticles.end(); it++ ) {
    out.push_back( &it->second );
  }
  std::sort( out.begin(), out.end(), compareParticlePtrs );
  return out;
}

Particles
Event::GetConstParticles() const
{
  Particles out;
  for ( ParticlesMap::const_iterator it=fParticles.begin(); it!=fParticles.end(); it++ ) {
    out.push_back( static_cast<Particle>( it->second ) );
  }
  std::sort( out.begin(), out.end(), compareParticle );
  return out;
}

ConstParticlesRef
Event::GetConstParticlesRef() const
{
  ConstParticlesRef out;
  for ( ParticlesMap::const_iterator it=fParticles.begin(); it!=fParticles.end(); it++ ) {
    out.push_back( &it->second );
  }
  std::sort( out.begin(), out.end(), compareParticlePtrs );
  return out;
}

ParticlesRef
Event::GetStableParticles()
{
  ParticlesRef out;
  for ( ParticlesMap::iterator it=fParticles.begin(); it!=fParticles.end(); it++ ) {
    if ( it->second.status==Particle::Undefined
      or it->second.status==Particle::FinalState ) {
      out.push_back( &it->second );
    }
  }
  std::sort( out.begin(), out.end() );
  return out;  
}

void
Event::Dump(bool stable_) const
{
  double pxtot, pytot, pztot, etot;
  std::ostringstream os;

  pxtot = pytot = pztot = etot = 0.;
  const ConstParticlesRef particles = GetConstParticlesRef();
  for ( ConstParticlesRef::const_iterator part_ref=particles.begin(); part_ref!=particles.end(); part_ref++ ) {
    const Particle* p = ( *part_ref );
    if ( stable_ and p->status!=Particle::FinalState ) continue;

    os << Form( "\n %2d\t%+6d", p->id, p->GetIntPDGId() );
    if ( p->name!="" ) os << Form( "%6s", p->name.c_str() );
    //else               os << std::setw(6) << Particle::ParticleCode( abs( p->GetPDGId() ) );
    else               os << "\t";
    os << "\t";
    if ( p->charge!=999. ) os << Form( "%6.2f\t", p->charge );
    else                   os << "\t";
    os << Form( "%4d\t%6d\t", p->role, p->status );
    if ( p->GetMothersIds().size()>0 ) 
      os << Form( "%2d(%2d)", *( p->GetMothersIds().begin() ), GetConstById( *( p->GetMothersIds().begin() ) ).role );
    else
      os << "      ";
    os << Form( "% 9.3f % 9.3f % 9.3f % 9.3f", p->GetMomentum().Px(), p->GetMomentum().Py(), p->GetMomentum().Pz(), p->E() );
    if ( p->status==Particle::Undefined
      or p->status==Particle::FinalState
      or p->status==Particle::Undecayed ) {
      const int sign = ( p->status==Particle::Undefined ) ? -1 : 1;
      pxtot += sign*p->GetMomentum().Px();
      pytot += sign*p->GetMomentum().Py();
      pztot += sign*p->GetMomentum().Pz();
      etot += sign*p->E();
    }
  }
  // We set a threshold to the computation precision
  if ( fabs(pxtot)<1.e-12 ) pxtot = 0.;
  if ( fabs(pytot)<1.e-12 ) pytot = 0.;
  if ( fabs(pztot)<1.e-12 ) pztot = 0.;
  if ( fabs( etot)<1.e-12 ) etot = 0.;
  //
  Information( Form( "Dump of event content:\n"
  "Part.\tPDG id\t\tCharge\tRole\tStatus\tMother\t\t4-Momentum (GeV)\n"
  "----\t------\t\t------\t----\t------\t------\t-------------------------------------"
  "%s\n"
  "---------------------------------------------------------------------------------------------\n"
  "Total:\t\t\t\t\t\t      % 9.4f % 9.4f % 9.4f % 9.4f", os.str().c_str(), pxtot, pytot, pztot, etot ) );
}
