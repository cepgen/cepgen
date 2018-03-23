#include "GenericProcess.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/ParticleProperties.h"

namespace CepGen
{
  namespace Process
  {
    const double GenericProcess::mp_ = ParticleProperties::mass( Proton );
    const double GenericProcess::mp2_ = GenericProcess::mp_*GenericProcess::mp_;

    GenericProcess::GenericProcess( const std::string& name, const std::string& description, bool has_event ) :
      first_run( true ),
      s_( 0. ), sqs_( 0. ), w1_( 0. ), w2_( 0. ), t1_( 0. ), t2_( 0. ), MX_( 0. ), MY_( 0. ),
      event_( new Event ),
      is_point_set_( false ), is_incoming_state_set_( false ), is_outgoing_state_set_( false ), is_kinematics_set_( false ),
      name_( name ), description_( description ), has_event_( has_event )
    {}

    GenericProcess::GenericProcess( const GenericProcess& proc ) :
      first_run( proc.first_run ),
      s_( proc.s_ ), sqs_( proc.sqs_ ),
      w1_( proc.w1_ ), w2_( proc.w2_ ),
      t1_( proc.w1_ ), t2_( proc.w2_ ),
      MX_( proc.w1_ ), MY_( proc.w2_ ),
      event_( new Event( *proc.event_.get() ) ),
      is_point_set_( proc.is_point_set_ ),
      is_incoming_state_set_( proc.is_incoming_state_set_ ), is_outgoing_state_set_( proc.is_outgoing_state_set_ ),
      is_kinematics_set_( proc.is_kinematics_set_ ),
      name_( proc.name_ ), description_( proc.description_ ),
      has_event_( proc.has_event_ )
    {}

    void
    GenericProcess::setPoint( const unsigned int ndim, double* x )
    {
      x_ = std::vector<double>( x, x+ndim );
      is_point_set_ = true;

      if ( Logger::get().level >= Logger::DebugInsideLoop )
        dumpPoint();
    }

    double
    GenericProcess::x( unsigned int idx ) const
    {
      if ( idx >= x_.size() )
        return -1.;
      return x_[idx];
    }

    void
    GenericProcess::prepareKinematics()
    {
      if ( !isKinematicsDefined() )
        throw Exception( __PRETTY_FUNCTION__, "Kinematics not properly defined for the process", FatalError );

      const Particle& ib1 = event_->getOneByRole( Particle::IncomingBeam1 );
      const Particle& ib2 = event_->getOneByRole( Particle::IncomingBeam2 );

      sqs_ = CMEnergy( ib1, ib2 );
      s_ = sqs_*sqs_;

      w1_ = ib1.mass2();
      w2_ = ib2.mass2();

      Debugging( Form( "Kinematics successfully prepared! sqrt(s) = %.2f", sqs_ ) );
    }

    void
    GenericProcess::dumpPoint() const
    {
      std::ostringstream os;
      for ( unsigned short i = 0; i < x_.size(); ++i ) {
        os << Form( "  x(%2d) = %8.6f\n\t", i, x_[i] );
      }
      Information( Form( "Number of integration parameters: %d\n\t"
                         "%s", x_.size(), os.str().c_str() ) );
    }

    void
    GenericProcess::setEventContent( const IncomingState& ini, const OutgoingState& fin )
    {
      event_->clear();
      //----- add the particles in the event

      //--- incoming state
      for ( const auto& ip : ini ) {
        Particle& p = event_->addParticle( ip.first );
        p.setPdgId( ip.second, ParticleProperties::charge( ip.second ) );
        p.setMass( ParticleProperties::mass( ip.second ) );
        if ( ip.first == Particle::IncomingBeam1
          || ip.first == Particle::IncomingBeam2 )
          p.setStatus( Particle::PrimordialIncoming );
      }
      //--- central system (if not already there)
      const auto& central_system = ini.find( Particle::CentralSystem );
      if ( central_system == ini.end() ) {
        Particle& p = event_->addParticle( Particle::Intermediate );
        p.setPdgId( invalidParticle );
        p.setStatus( Particle::Propagator );
      }
      //--- outgoing state
      for ( const auto& opl : fin ) { // pair(role, list of PDGids)
        for ( const auto& pdg : opl.second ) {
          Particle& p = event_->addParticle( opl.first );
          p.setPdgId( pdg, ParticleProperties::charge( pdg ) );
          p.setMass( ParticleProperties::mass( pdg ) );
        }
      }

      //----- define the particles parentage

      const Particles parts = event_->particles();
      for ( const auto& p : parts ) {
        Particle& part = event_->getById( p.id() );
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

    Particles&
    GenericProcess::particles( const Particle::Role& role )
    {
      return event_->getByRole( role );
    }

    bool
    GenericProcess::isKinematicsDefined()
    {
      // check the incoming state
      if ( !particles( Particle::IncomingBeam1 ).empty()
        && !particles( Particle::IncomingBeam2 ).empty() )
        is_incoming_state_set_ = true;

      // check the outgoing state
      if ( !particles( Particle::OutgoingBeam1 ).empty()
        && !particles( Particle::OutgoingBeam2 ).empty()
        && !particles( Particle::CentralSystem ).empty() )
        is_outgoing_state_set_ = true;

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
