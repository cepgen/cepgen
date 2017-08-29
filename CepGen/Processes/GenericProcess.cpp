#include "GenericProcess.h"

namespace CepGen
{
  namespace Process
  {
    GenericProcess::GenericProcess( const std::string& name, const std::string& description, bool has_event ) :
      event_( std::shared_ptr<Event>( new Event ) ),
      is_point_set_( false ), is_incoming_state_set_( false ), is_outgoing_state_set_( false ), is_kinematics_set_( false ),
      name_( name ), description_( description ),
      total_gen_time_( 0. ), num_gen_events_( 0 ), has_event_( has_event )
    {}

    GenericProcess::~GenericProcess()
    {}

    void
    GenericProcess::setPoint( const unsigned int ndim, double* x )
    {
      if ( ndim != x_.size() ) x_.resize( ndim );

      x_ = std::vector<double>( x, x+ndim );
      is_point_set_ = true;
      if ( Logger::get().level>=Logger::DebugInsideLoop ) { dumpPoint( DebugMessage ); }
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
    GenericProcess::dumpPoint( const ExceptionType& et )
    {
    std::ostringstream os;
    for ( unsigned int i=0; i<x_.size(); i++ ) {
      os << Form( "  x(%2d) = %8.6f\n\t", i, x_[i] );
    }
    if ( et < DebugMessage ) {
      Information( Form( "Number of integration parameters: %d\n\t"
                         "%s", x_.size(), os.str().c_str() ) ); }
    else {
      Debugging( Form( "Number of integration parameters: %d\n\t"
                       "%s", x_.size(), os.str().c_str() ) ); }
    }

    void
    GenericProcess::setEventContent( const IncomingState& is, const OutgoingState& os )
    {  
      for ( IncomingState::const_iterator ip=is.begin(); ip!=is.end(); ip++ ) {
        event_->addParticle( Particle( ip->first, ip->second ), true );
      }

      // Prepare the central system if not already there
      IncomingState::const_iterator central_system = is.find( Particle::CentralSystem );
      if ( central_system==is.end() ) {
        Particle cs( Particle::Intermediate, Particle::invalidParticle );
        cs.addMother( event_->getOneByRole( Particle::Parton1 ) );
        cs.addMother( event_->getOneByRole( Particle::Parton2 ) );
        cs.setStatus( Particle::Propagator );
        event_->addParticle( cs, true );
      }

      for ( OutgoingState::const_iterator op=os.begin(); op!=os.end(); op++ ) {
        for ( std::vector<Particle::ParticleCode>::const_iterator it=op->second.begin(); it!=op->second.end(); ++it ) {
          event_->addParticle( Particle( op->first, *it ), false );
        }
      }
  
      // Incoming particles (incl. eventual partons)
      for ( IncomingState::const_iterator ip=is.begin(); ip!=is.end(); ip++ ) {
        Particle& p = event_->getOneByRole( ip->first );
        p.setStatus( Particle::Undefined );
        switch ( ip->first ) {
          case Particle::IncomingBeam1:
          case Particle::IncomingBeam2: break;
          case Particle::Parton1:      p.addMother( event_->getOneByRole( Particle::IncomingBeam1 ) ); break;
          case Particle::Parton2:      p.addMother( event_->getOneByRole( Particle::IncomingBeam2 ) ); break;
          case Particle::Intermediate: p.addMother( event_->getOneByRole( Particle::Parton1 ) ); break;
          default: break;
        }
      }
      // Outgoing particles (central, and outgoing primary particles or remnants)
      for ( OutgoingState::const_iterator op=os.begin(); op!=os.end(); op++ ) {
        Particles& parts = event_->getByRole( op->first );
        for ( Particles::iterator p = parts.begin(); p != parts.end(); ++p ) {
          p->setStatus( Particle::Undefined );
          switch ( op->first ) {
            case Particle::OutgoingBeam1: p->addMother( event_->getOneByRole( Particle::IncomingBeam1 ) ); break;
            case Particle::OutgoingBeam2: p->addMother( event_->getOneByRole( Particle::IncomingBeam2 ) ); break;
            case Particle::CentralSystem: p->addMother( event_->getOneByRole( Particle::Intermediate ) ); break;
            default: break;
          }
        }
      }
      event_->init();
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

      bool inel_p1 = false,
           inel_p2 = false;

      switch ( cuts_.mode ) {
        case Kinematics::ElectronElectron: {
          fp1 = TrivialFormFactors(); // electron (trivial) form factor
          fp2 = TrivialFormFactors(); // electron (trivial) form factor
        } break;
        case Kinematics::ProtonElectron: {
          fp1 = ElasticFormFactors( -t1_, w1_ ); // proton elastic form factor
          fp2 = TrivialFormFactors(); // electron (trivial) form factor
        } break;
        case Kinematics::ElectronProton: {
          fp1 = TrivialFormFactors(); // electron (trivial) form factor
          fp2 = ElasticFormFactors( -t2_, w2_ ); // proton elastic form factor
        } break;
        case Kinematics::ElasticElastic: {
          fp1 = ElasticFormFactors( -t1_, w1_ ); // proton elastic form factor
          fp2 = ElasticFormFactors( -t2_, w2_ ); // proton elastic form factor
        } break;
        case Kinematics::ElasticInelastic: {
          fp1 = ElasticFormFactors( -t1_, w1_ );
          inel_p2 = true;
        } break;
        case Kinematics::InelasticElastic: {
          inel_p1 = true;
          fp2 = ElasticFormFactors( -t2_, w2_ );
        } break;
        case Kinematics::InelasticInelastic: {
          inel_p1 = inel_p2 = true;
        } break;
      }
      switch ( cuts_.remnant_mode ) {
        case SuriYennie:
        default: {
          if ( inel_p1 ) fp1 = SuriYennieFormFactors( -t1_, w1_, mx2 );
          if ( inel_p2 ) fp2 = SuriYennieFormFactors( -t2_, w2_, my2 );
        } break;
        case Fiore:
        case FioreSea:
        case FioreVal: { // low-Q2 inelastic form factor
          if ( inel_p1 ) fp1 = FioreBrasseFormFactors( -t1_, w1_, mx2 );
          if ( inel_p2 ) fp2 = FioreBrasseFormFactors( -t2_, w2_, my2 );
        } break;
        case SzczurekUleshchenko: {
          if ( inel_p1 ) fp1 = SzczurekUleshchenkoFormFactors( -t1_, w1_, mx2 );
          if ( inel_p2 ) fp2 = SzczurekUleshchenkoFormFactors( -t2_, w2_, my2 );
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
