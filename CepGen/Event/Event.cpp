#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include <algorithm>
#include <math.h>

namespace cepgen
{
  Event::Event( bool compressed ) :
    num_hadronisation_trials( 0 ),
    time_generation( -1. ), time_total( -1. ), weight( 0. ),
    compressed_( compressed )
  {}

  Event::Event( const Event& rhs ) :
    num_hadronisation_trials( rhs.num_hadronisation_trials ),
    time_generation( rhs.time_generation ), time_total( rhs.time_total ),
    weight( rhs.weight ),
    particles_( rhs.particles_ ),
    evtcontent_( rhs.evtcontent_ ),
    compressed_( rhs.compressed_ )
  {}

  void
  Event::clear()
  {
    particles_.clear();
    time_generation = -1.;
    time_total = -1.;
    weight = 0.;
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

  bool
  Event::compressed() const
  {
    return compressed_;
  }

  Event
  Event::compress() const
  {
    if ( compressed_ )
      return *this;
    Event out( /*compressed=*/true );
    size_t i = 0;
    //--- add all necessary particles
    for ( const auto& role : { Particle::IncomingBeam1, Particle::IncomingBeam2,
                               Particle::OutgoingBeam1, Particle::OutgoingBeam2,
                               Particle::Parton1, Particle::Parton2,
                               Particle::CentralSystem } ) {
      for ( const auto& old_part : operator[]( role ) ) {
        auto& new_part = out.addParticle( role );
        new_part = old_part; // copy all attributes
        new_part.setId( i++ );
        new_part.clearMothers();
        new_part.clearDaughters();
      }
    }
    //--- fix parentage for outgoing beam particles
    if ( out[Particle::OutgoingBeam1].size() > 1
      || out[Particle::OutgoingBeam2].size() > 1 )
      CG_WARNING( "Event:compress" )
        << "Event compression not designed for already fragmented beam remnants!\n\t"
        << "Particles parentage is not guaranteed to be conserved.";
    for ( auto& part : out[Particle::OutgoingBeam1] )
      part.addMother( out[Particle::IncomingBeam1][0] );
    for ( auto& part : out[Particle::OutgoingBeam2] )
      part.addMother( out[Particle::IncomingBeam2][0] );
    //--- fix parentage for incoming partons
    for ( auto& part : out[Particle::Parton1] )
      part.addMother( out[Particle::IncomingBeam1][0] );
    for ( auto& part : out[Particle::Parton2] )
      part.addMother( out[Particle::IncomingBeam2][0] );
    //--- fix parentage for central system
    for ( auto& part : out[Particle::CentralSystem] ) {
      part.addMother( out[Particle::Parton1][0] );
      part.addMother( out[Particle::Parton2][0] );
    }
    return out;
  }

  double
  Event::cmEnergy() const
  {
    return CMEnergy( oneWithRole( Particle::IncomingBeam1 ), oneWithRole( Particle::IncomingBeam2 ) );
  }

  Particles&
  Event::operator[]( Particle::Role role )
  {
    //--- retrieve all particles with a given role
    return particles_[role];
  }

  const Particles&
  Event::operator[]( Particle::Role role ) const
  {
    if ( particles_.count( role ) == 0 )
      throw CG_FATAL( "Event" ) << "Failed to retrieve a particle with " << role << " role.";
    //--- retrieve all particles with a given role
    return particles_.at( role );
  }

  ParticlesIds
  Event::ids( Particle::Role role ) const
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
  Event::oneWithRole( Particle::Role role )
  {
    //--- retrieve the first particle of a given role
    Particles& parts_by_role = operator[]( role );
    if ( parts_by_role.empty() )
      throw CG_FATAL( "Event" ) << "No particle retrieved with " << role << " role.";
    if ( parts_by_role.size() > 1 )
      throw CG_FATAL( "Event" ) << "More than one particle with " << role << " role: "
        << parts_by_role.size() << " particles.";
    return *parts_by_role.begin();
  }

  const Particle&
  Event::oneWithRole( Particle::Role role ) const
  {
    //--- retrieve the first particle of a given role
    const Particles& parts_by_role = operator[]( role );
    if ( parts_by_role.empty() )
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
  Event::operator[]( int id ) const
  {
    for ( const auto& role_part : particles_ )
      for ( const auto& part : role_part.second )
        if ( part.id() == id )
          return part;

    throw CG_FATAL( "Event" ) << "Failed to retrieve the particle with id=" << id << ".";
  }

  Particles
  Event::operator[]( const ParticlesIds& ids ) const
  {
    Particles out;
    for ( const auto& id : ids )
      out.emplace_back( operator[]( id ) );

    return out;
  }

  Particles
  Event::mothers( const Particle& part ) const
  {
    return operator[]( part.mothers() );
  }

  Particles
  Event::daughters( const Particle& part ) const
  {
    return operator[]( part.daughters() );
  }

  ParticleRoles
  Event::roles() const
  {
    ParticleRoles out;
    for ( const auto& pr : particles_ )
      out.emplace_back( pr.first );
    return out;
  }

  Particle&
  Event::addParticle( Particle& part, bool replace )
  {
    CG_DEBUG_LOOP( "Event" ) << "Particle with PDGid = " << part.integerPdgId() << " has role " << part.role();
    if ( part.role() <= 0 )
      throw CG_FATAL( "Event" ) << "Trying to add a particle with role=" << (int)part.role() << ".";

    //--- retrieve the list of particles with the same role
    Particles& part_with_same_role = operator[]( part.role() );

    //--- specify the id
    if ( part_with_same_role.empty() && part.id() < 0 ) part.setId( size() ); // set the id if previously invalid/inexistent
    if ( !part_with_same_role.empty() ) {
      if ( replace ) part.setId( part_with_same_role[0].id() ); // set the previous id if replacing a particle
      else part.setId( size() );
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
  Event::size() const
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
      Momentum ptot;
      for ( const auto& daugh : daughters ) {
        const Particle& d = operator[]( daugh );
        const ParticlesIds mothers = d.mothers();
        ptot += d.momentum();
        if ( mothers.size() < 2 )
          continue;
        for ( const auto& moth : mothers )
          if ( moth != part.id() )
            ptot -= operator[]( moth ).momentum();
      }
      const double mass_diff = ( ptot-part.momentum() ).mass();
      if ( fabs( mass_diff ) > MIN_PRECISION ) {
        dump();
        throw CG_FATAL( "Event" ) << "Error in momentum balance for particle " << part.id() << ": mdiff = " << mass_diff << ".";
      }
    }
  }

  void
  Event::dump() const
  {
    CG_INFO( "Event" ) << *this;
  }

  std::ostream&
  operator<<( std::ostream& out, const Event& ev )
  {
    const Particles parts = ev.particles();
    std::ostringstream os;

    Momentum p_total;
    for ( const auto& part : parts ) {
      const ParticlesIds mothers = part.mothers();
      {
        std::ostringstream oss_pdg;
        if ( part.pdgId() == PDG::invalid && !mothers.empty() ) {
          //--- if particles compound
          std::string delim;
          for ( unsigned short i = 0; i < mothers.size(); ++i )
            try {
              oss_pdg << delim
                << PDG::get().name( ev[*std::next( mothers.begin(), i )].pdgId() ), delim = "/";
            } catch ( const Exception& ) {
              oss_pdg << delim
                << ev[*std::next( mothers.begin(), i )].pdgId(), delim = "/";
            }
          os << utils::format( "\n %2d\t\t   %-7s", part.id(), oss_pdg.str().c_str() );
        }
        else {
          //--- if single particle/HI
          if ( (HeavyIon)part.pdgId() )
            oss_pdg << (HeavyIon)part.pdgId();
          else
            try {
              oss_pdg << PDG::get().name( part.pdgId() );
            } catch ( const Exception& ) {
              oss_pdg << "?";
            }
          os << utils::format( "\n %2d\t%-+10d %-7s", part.id(), part.integerPdgId(), oss_pdg.str().c_str() );
        }
      }
      os << "\t";
      if ( part.charge() != (int)part.charge() ) {
        if ( part.charge()*2 == (int)( part.charge()*2 ) )
          os << utils::format( "%-d/2", (int)( part.charge()*2 ) );
        else if ( part.charge()*3 == (int)( part.charge()*3 ) )
          os << utils::format( "%-d/3", (int)( part.charge()*3 ) );
        else
          os << utils::format( "%-.2f", part.charge() );
      }
      else
        os << utils::format( "%-g", part.charge() );
      os << "\t";
      { std::ostringstream oss; oss << part.role(); os << utils::format( "%-8s %6d\t", oss.str().c_str(), part.status() ); }
      if ( !mothers.empty() ) {
        std::ostringstream oss;
        unsigned short i = 0;
        for ( const auto& moth : mothers ) {
          oss << ( i > 0 ? "+" : "" ) << moth;
          ++i;
        }
        os << utils::format( "%6s ", oss.str().c_str() );
      }
      else os << "       ";
      const auto& mom = part.momentum();
      os << utils::format( "% 9.6e % 9.6e % 9.6e % 9.6e % 12.5f", mom.px(), mom.py(), mom.pz(), part.energy(), part.mass() );

      // discard non-primary, decayed particles
      if ( part.status() >= Particle::Status::Undefined ) {
        const int sign = ( part.status() == Particle::Status::Undefined )
          ? -1
          : +1;
        p_total += sign*mom;
      }
    }
    //--- set a threshold to the computation precision
    p_total.truncate();
    //
    return out << utils::format(
       "Event content:\n"
       " Id\tPDG id\t   Name\t\tCharge\tRole\t Status\tMother\tpx            py            pz            E      \t M         \n"
       " --\t------\t   ----\t\t------\t----\t ------\t------\t----GeV/c---  ----GeV/c---  ----GeV/c---  ----GeV/c---\t --GeV/c²--"
       "%s\n"
       " ----------------------------------------------------------------------------------------------------------------------------------\n"
       "\t\t\t\t\t\t\tBalance% 9.6e % 9.6e % 9.6e % 9.6e", os.str().c_str(), p_total.px(), p_total.py(), p_total.pz(), p_total.energy() );
  }

  //------------------------------------------------------------------------------------------------

  Event::NumParticles::NumParticles() :
    cs( 0 ), op1( 0 ), op2( 0 )
  {}

  Event::NumParticles::NumParticles( const NumParticles& np ) :
    cs( np.cs ), op1( np.op1 ), op2( np.op2 )
  {}
}

