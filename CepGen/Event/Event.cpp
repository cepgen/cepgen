#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <algorithm>
#include <math.h>

namespace CepGen
{
  Event::Event() :
    num_hadronisation_trials( 0 ),
    time_generation( -1. ), time_total( -1. )
  {}

  Event::Event( const Event& rhs ) :
    num_hadronisation_trials( rhs.num_hadronisation_trials ),
    time_generation( rhs.time_generation ), time_total( rhs.time_total ),
    particles_( rhs.particles_ ),
    evtcontent_( rhs.evtcontent_ )
  {}

  void
  Event::clear()
  {
    particles_.clear();
    time_generation = -1.;
    time_total = -1.;
  }

  void
  Event::freeze()
  {
    //--- store a snapshot of the primordial event block
    if ( particles_.count( Particle::CentralSystem ) > 0 )
      evtcontent_.cs = particles_[Particle::CentralSystem].size();
    if ( particles_.count( Particle::OutgoingBeam1 ) > 0 )
      evtcontent_.op1 = particles_[Particle::OutgoingBeam1].size();
    if ( particles_.count( Particle::OutgoingBeam2 ) > 0 )
      evtcontent_.op2 = particles_[Particle::OutgoingBeam2].size();
  }

  void
  Event::restore()
  {
    //--- remove all particles after the primordial event block
    if ( particles_.count( Particle::CentralSystem ) > 0 )
      particles_[Particle::CentralSystem].resize( evtcontent_.cs );
    if ( particles_.count( Particle::OutgoingBeam1 ) > 0 )
      particles_[Particle::OutgoingBeam1].resize( evtcontent_.op1 );
    if ( particles_.count( Particle::OutgoingBeam2 ) > 0 )
      particles_[Particle::OutgoingBeam2].resize( evtcontent_.op2 );
  }

  double
  Event::cmEnergy() const
  {
    return CMEnergy( getOneByRole( Particle::IncomingBeam1 ), getOneByRole( Particle::IncomingBeam2 ) );
  }

  Particles&
  Event::getByRole( Particle::Role role )
  {
    //--- retrieve all particles with a given role
    return particles_[role];
  }

  const Particles&
  Event::getByRole( Particle::Role role ) const
  {
    if ( particles_.count( role ) == 0 )
      throw CG_FATAL( "Event" ) << "Failed to retrieve a particle with " << role << " role.";
    //--- retrieve all particles with a given role
    return particles_.at( role );
  }

  ParticlesIds
  Event::getIdsByRole( Particle::Role role ) const
  {
    ParticlesIds out;
    //--- retrieve all particles ids with a given role
    if ( particles_.count( role ) == 0 )
      return out;

    for ( const auto& part : particles_.at( role ) )
      out.insert( part.id() );

    return out;
  }

  Particle&
  Event::getOneByRole( Particle::Role role )
  {
    //--- retrieve the first particle a the given role
    Particles& parts_by_role = getByRole( role );
    if ( parts_by_role.size() == 0 )
      throw CG_FATAL( "Event" ) << "No particle retrieved with " << role << " role.";
    if ( parts_by_role.size() > 1 )
      throw CG_FATAL( "Event" ) << "More than one particle with " << role << " role: "
        << parts_by_role.size() << " particles.";
    return *parts_by_role.begin();
  }

  const Particle&
  Event::getOneByRole( Particle::Role role ) const
  {
    if ( particles_.count( role ) == 0 )
      throw CG_FATAL( "Event" ) << "Failed to retrieve a particle with " << role << " role.";
    //--- retrieve the first particle a the given role
    const Particles& parts_by_role = particles_.at( role );
    if ( parts_by_role.size() == 0 )
      throw CG_FATAL( "Event" ) << "No particle retrieved with " << role << " role.";
    if ( parts_by_role.size() > 1 )
      throw CG_FATAL( "Event" ) << "More than one particle with " << role << " role: "
        << parts_by_role.size() << " particles";
    return *parts_by_role.begin();
  }

  Particle&
  Event::operator[]( int id )
  {
    for ( auto& role_part : particles_ )
      for ( auto& part : role_part.second )
        if ( part.id() == id )
          return part;

    throw CG_FATAL( "Event" ) << "Failed to retrieve the particle with id=" << id << ".";
  }

  const Particle&
  Event::getConstById( int id ) const
  {
    for ( const auto& role_part : particles_ )
      for ( const auto& part : role_part.second )
        if ( part.id() == id )
          return part;

    throw CG_FATAL( "Event" ) << "Failed to retrieve the particle with id=" << id << ".";
  }

  Particles
  Event::getByIds( const ParticlesIds& ids ) const
  {
    Particles out;
    for ( const auto& id : ids )
      out.emplace_back( getConstById( id ) );

    return out;
  }

  Particles
  Event::mothers( const Particle& part )
  {
    return getByIds( part.mothers() );
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
    ParticlesMap::const_iterator it, end = particles_.end();
    for ( it = particles_.begin(); it != end; it = particles_.upper_bound( it->first ) ) {
      out.emplace_back( it->first );
    }
    return out;
  }

  Particle&
  Event::addParticle( Particle& part, bool replace )
  {
    CG_DEBUG_LOOP( "Event" ) << "Particle with PDGid = " << part.integerPdgId() << " has role " << part.role();
    if ( part.role() <= 0 )
      throw CG_FATAL( "Event" ) << "Trying to add a particle with role=" << (int)part.role() << ".";

    //--- retrieve the list of particles with the same role
    Particles& part_with_same_role = getByRole( part.role() );

    //--- specify the id
    if ( part_with_same_role.empty() && part.id() < 0 ) part.setId( numParticles() ); // set the id if previously invalid/inexistent
    if ( !part_with_same_role.empty() ) {
      if ( replace ) part.setId( part_with_same_role[0].id() ); // set the previous id if replacing a particle
      else part.setId( numParticles() );
    }

    //--- add the particle to the collection
    if ( replace ) part_with_same_role = Particles( 1, part ); // generate a vector containing only this particle
    else part_with_same_role.emplace_back( part );

    return part_with_same_role.back();
  }

  Particle&
  Event::addParticle( Particle::Role role, bool replace )
  {
    Particle np( role, PDG::invalid );
    return addParticle( np, replace );
  }

  size_t
  Event::numParticles() const
  {
    size_t out = 0;
    for ( const auto& role_part : particles_ )
      out += role_part.second.size();
    return out;
  }

  const Particles
  Event::particles() const
  {
    Particles out;
    for ( const auto& role_part : particles_ )
      out.insert( out.end(), role_part.second.begin(), role_part.second.end() );

    std::sort( out.begin(), out.end() );
    return out;
  }

  const Particles
  Event::stableParticles() const
  {
    Particles out;
    for ( const auto& role_part : particles_ )
      for ( const auto& part : role_part.second )
        if ( (short)part.status() > 0 )
          out.emplace_back( part );

    std::sort( out.begin(), out.end() );
    return out;
  }

  void
  Event::checkKinematics() const
  {
    // check the kinematics through parentage
    for ( const auto& part : particles() ) {
      ParticlesIds daughters = part.daughters();
      if ( daughters.empty() )
        continue;
      Particle::Momentum ptot;
      for ( const auto& daugh : daughters ) {
        const Particle& d = getConstById( daugh );
        const ParticlesIds mothers = d.mothers();
        ptot += d.momentum();
        if ( mothers.size() < 2 )
          continue;
        for ( const auto& moth : mothers ) {
          if ( moth == part.id() )
            continue;
          ptot -= getConstById( moth ).momentum();
        }
      }
      const double mass_diff = ( ptot-part.momentum() ).mass();
      if ( fabs( mass_diff ) > minimal_precision_ ) {
        dump();
        throw CG_FATAL( "Event" ) << "Error in momentum balance for particle " << part.id() << ": mdiff = " << mass_diff << ".";
      }
    }
  }

  void
  Event::dump( std::ostream& out, bool stable ) const
  {
    const Particles parts = ( stable ) ? stableParticles() : particles();

    std::ostringstream os;

    double pxtot = 0., pytot = 0., pztot = 0., etot = 0.;
    for ( const auto& part : parts ) {
      const ParticlesIds mothers = part.mothers();
      {
        std::ostringstream oss_pdg;
        if ( part.pdgId() == PDG::invalid && mothers.size() > 0 ) {
          for ( unsigned short i = 0; i < mothers.size(); ++i )
            oss_pdg << ( i > 0 ? "/" : "" ) << getConstById( *std::next( mothers.begin(), i ) ).pdgId();
          os << Form( "\n %2d\t\t%-10s", part.id(), oss_pdg.str().c_str() );
        }
        else {
          if ( (HeavyIon)part.pdgId() )
            oss_pdg << (HeavyIon)part.pdgId();
          else
            oss_pdg << part.pdgId();
          os << Form( "\n %2d\t%-+7d %-10s", part.id(), part.integerPdgId(), oss_pdg.str().c_str() );
        }
      }
      os << "\t";
      if ( part.charge() != 999. )
        os << Form( "%-.2f\t", part.charge() );
      else
        os << "\t";
      { std::ostringstream oss; oss << part.role(); os << Form( "%-8s %6d\t", oss.str().c_str(), part.status() ); }
      if ( !mothers.empty() ) {
        std::ostringstream oss;
        unsigned short i = 0;
        for ( const auto& moth : mothers ) {
          oss << ( i > 0 ? "+" : "" ) << moth;
          ++i;
        }
        os << Form( "%6s ", oss.str().c_str() );
      }
      else os << "       ";
      const Particle::Momentum mom = part.momentum();
      os << Form( "% 9.6e % 9.6e % 9.6e % 9.6e % 12.5f", mom.px(), mom.py(), mom.pz(), part.energy(), part.mass() );

      // discard non-primary, decayed particles
      if ( (short)part.status() >= 0. ) {
        const int sign = ( part.status() == Particle::Status::Undefined ) ? -1 : 1;
        pxtot += sign*mom.px();
        pytot += sign*mom.py();
        pztot += sign*mom.pz();
        etot += sign*part.energy();
      }
    }
    //--- set a threshold to the computation precision
    if ( fabs( pxtot ) < minimal_precision_ ) pxtot = 0.;
    if ( fabs( pytot ) < minimal_precision_ ) pytot = 0.;
    if ( fabs( pztot ) < minimal_precision_ ) pztot = 0.;
    if ( fabs(  etot ) < minimal_precision_ ) etot = 0.;
    //
    CG_INFO( "Event" )
     << Form( "Dump of event content:\n"
              " Id\tPDG id\tName\t\tCharge\tRole\t Status\tMother\tpx            py            pz            E      \t M         \n"
              " --\t------\t----\t\t------\t----\t ------\t------\t----GeV/c---  ----GeV/c---  ----GeV/c---  ----GeV/c---\t --GeV/c²--"
              "%s\n"
              " ----------------------------------------------------------------------------------------------------------------------------------\n"
              "\t\t\t\t\t\t\tBalance% 9.6e % 9.6e % 9.6e % 9.6e", os.str().c_str(), pxtot, pytot, pztot, etot );
  }

  //------------------------------------------------------------------------------------------------

  Event::NumParticles::NumParticles() :
    cs( 0 ), op1( 0 ), op2( 0 )
  {}

  Event::NumParticles::NumParticles( const NumParticles& np ) :
    cs( np.cs ), op1( np.op1 ), op2( np.op2 )
  {}
}
