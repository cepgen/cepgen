#include "GenericProcess.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace CepGen
{
  namespace Process
  {
    GenericProcess::GenericProcess( const std::string& name, const std::string& description, bool has_event ) :
      s_( 0. ), sqs_( 0. ), w1_( 0. ), w2_( 0. ), t1_( 0. ), t2_( 0. ), MX_( 0. ), MY_( 0. ),
      event_( std::shared_ptr<Event>( new Event ) ),
      is_point_set_( false ), is_incoming_state_set_( false ), is_outgoing_state_set_( false ), is_kinematics_set_( false ),
      name_( name ), description_( description ),
      total_gen_time_( 0. ), num_gen_events_( 0 ), has_event_( has_event )
    {}

    void
    GenericProcess::setPoint( const unsigned int ndim, double* x )
    {
      if ( ndim != x_.size() ) x_.resize( ndim );

      x_ = std::vector<double>( x, x+ndim );
      is_point_set_ = true;
      if ( Logger::get().level>=Logger::DebugInsideLoop ) { dumpPoint(); }
    }

    void
    GenericProcess::prepareKinematics()
    {
      if ( !isKinematicsDefined() ) return; // FIXME dump some information...
      const Particle ib1 = event_->getOneByRole( Particle::IncomingBeam1 ),
                     ib2 = event_->getOneByRole( Particle::IncomingBeam2 );

      sqs_ = CMEnergy( ib1, ib2 );
      s_ = sqs_*sqs_;

      w1_ = ib1.mass2();
      w2_ = ib2.mass2();

      Debugging( Form( "Kinematics successfully prepared! sqrt(s) = %.2f", sqs_ ) );
    }

    void
    GenericProcess::dumpPoint()
    {
      std::ostringstream os;
      for ( unsigned int i = 0; i < x_.size(); ++i ) {
        os << Form( "  x(%2d) = %8.6f\n\t", i, x_[i] );
      }
      Information( Form( "Number of integration parameters: %d\n\t"
                         "%s", x_.size(), os.str().c_str() ) );
    }

    void
    GenericProcess::setEventContent( const IncomingState& is, const OutgoingState& os )
    {
      event_->clear();
      //----- add the particles in the event

      //--- incoming state
      for ( IncomingState::const_iterator ip = is.begin(); ip != is.end(); ++ip ) {
        Particle& p = event_->addParticle( ip->first );
        p.setPdgId( ip->second, ParticleProperties::charge( ip->second ) );
      }
      //--- central system (if not already there)
      IncomingState::const_iterator central_system = is.find( Particle::CentralSystem );
      if ( central_system == is.end() ) {
        Particle& p = event_->addParticle( Particle::Intermediate );
        p.setPdgId( invalidParticle );
        p.setStatus( Particle::Propagator );
      }
      //--- outgoing state
      for ( OutgoingState::const_iterator op = os.begin(); op != os.end(); ++op ) {
        for ( std::vector<ParticleCode>::const_iterator it = op->second.begin(); it != op->second.end(); ++it ) {
          Particle& p = event_->addParticle( op->first );
          p.setPdgId( *it, ParticleProperties::charge( *it ) );
        }
      }

      //----- define the particles parentage

      const Particles parts = event_->particles();
      for ( Particles::const_iterator p = parts.begin(); p != parts.end(); ++p ) {
        Particle& part = event_->getById( p->id() );
        switch ( part.role() ) {
          case Particle::OutgoingBeam1:
          case Particle::Parton1:
            part.addMother( event_->getOneByRole( Particle::IncomingBeam1 ) );
            break;
          case Particle::OutgoingBeam2:
          case Particle::Parton2:
            part.addMother( event_->getOneByRole( Particle::IncomingBeam2 ) );
            break;
          case Particle::Intermediate:
            part.addMother( event_->getOneByRole( Particle::Parton1 ) );
            part.addMother( event_->getOneByRole( Particle::Parton2 ) );
            break;
          case Particle::CentralSystem:
            part.addMother( event_->getOneByRole( Particle::Intermediate ) );
            break;
          default: break;
        }
      }

      //----- freeze the event as it is

      event_->freeze();
    }

    void
    GenericProcess::setIncomingKinematics( const Particle::Momentum& p1, const Particle::Momentum& p2 )
    {
      event_->getOneByRole( Particle::IncomingBeam1 ).setMomentum( p1 );
      event_->getOneByRole( Particle::IncomingBeam2 ).setMomentum( p2 );
    }

    void
    GenericProcess::formFactors( double q1, double q2, FormFactors& fp1, FormFactors& fp2 ) const
    {
      const double mx2 = MX_*MX_, my2 = MY_*MY_;

      switch ( cuts_.mode ) {
        case Kinematics::ElectronElectron: {
          fp1 = FormFactors::Trivial(); // electron (trivial) form factor
          fp2 = FormFactors::Trivial(); // electron (trivial) form factor
        } break;
        case Kinematics::ProtonElectron: {
          fp1 = FormFactors::ProtonElastic( -t1_ ); // proton elastic form factor
          fp2 = FormFactors::Trivial(); // electron (trivial) form factor
        } break;
        case Kinematics::ElectronProton: {
          fp1 = FormFactors::Trivial(); // electron (trivial) form factor
          fp2 = FormFactors::ProtonElastic( -t2_ ); // proton elastic form factor
        } break;
        case Kinematics::ElasticElastic: {
          fp1 = FormFactors::ProtonElastic( -t1_ ); // proton elastic form factor
          fp2 = FormFactors::ProtonElastic( -t2_ ); // proton elastic form factor
        } break;
        case Kinematics::ElasticInelastic: {
          fp1 = FormFactors::ProtonElastic( -t1_ );
          fp2 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t2_, w2_, my2 );
        } break;
        case Kinematics::InelasticElastic: {
          fp1 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t1_, w1_, mx2 );
          fp2 = FormFactors::ProtonElastic( -t2_ );
        } break;
        case Kinematics::InelasticInelastic: {
          fp1 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t1_, w1_, mx2 );
          fp2 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t2_, w2_, my2 );
        } break;
      }
    }

    inline bool
    GenericProcess::isKinematicsDefined()
    {
      // check the incoming state
      if ( !particles( Particle::IncomingBeam1 ).empty() && !particles( Particle::IncomingBeam2 ).empty() ) {
        is_incoming_state_set_ = true;
      }
      // check the outgoing state
      if ( !particles( Particle::OutgoingBeam1 ).empty()
        && !particles( Particle::OutgoingBeam2 ).empty()
        && !particles( Particle::CentralSystem ).empty() ) {
        is_outgoing_state_set_ = true;
      }
      // combine both states
      is_kinematics_set_ = is_incoming_state_set_ && is_outgoing_state_set_;
      return is_kinematics_set_;
    }

    std::ostream&
    operator<<( std::ostream& os, const GenericProcess& proc )
    {
      return os << proc.name().c_str();
    }

    std::ostream&
    operator<<( std::ostream& os, const GenericProcess* proc )
    {
      return os << proc->name().c_str();
    }
  }
}
