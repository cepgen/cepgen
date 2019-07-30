#include "CepGen/Processes/GenericProcess.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen
{
  namespace proc
  {
    const double GenericProcess::mp_ = PDG::get().mass( PDG::proton );
    const double GenericProcess::mp2_ = GenericProcess::mp_*GenericProcess::mp_;

    GenericProcess::GenericProcess( const ParametersList& params, const std::string& name, const std::string& description, bool has_event ) :
      params_( params ), name_( name ), description_( description ),
      first_run( true ),
      s_( -1. ), sqs_( -1. ),
      MX_( -1. ), MY_( -1. ), w1_( -1. ), w2_( -1. ),
      t1_( -1. ), t2_( -1. ),
      mode_( (KinematicsMode)params.get<int>( "mode" ) ),
      has_event_( has_event ), event_( new Event ),
      is_point_set_( false )
    {}

    GenericProcess::GenericProcess( const GenericProcess& proc ) :
      params_( proc.params_ ), name_( proc.name_ ), description_( proc.description_ ),
      first_run( proc.first_run ),
      s_( proc.s_ ), sqs_( proc.sqs_ ),
      MX_( proc.MX_ ), MY_( proc.MY_ ), w1_( proc.w1_ ), w2_( proc.w2_ ),
      t1_( -1. ), t2_( -1. ),
      mode_( proc.mode_ ), kin_( proc.kin_ ),
      has_event_( proc.has_event_ ), event_( new Event( *proc.event_.get() ) ),
      is_point_set_( false )
    {}

    GenericProcess&
    GenericProcess::operator=( const GenericProcess& proc )
    {
      params_ = proc.params_;
      name_ = proc.name_; description_ = proc.description_;
      first_run = proc.first_run;
      s_ = proc.s_; sqs_ = proc.sqs_;
      MX_ = proc.MX_; MY_ = proc.MY_; w1_ = proc.w1_; w2_ = proc.w2_;
      mode_ = proc.mode_;
      kin_ = proc.kin_;
      has_event_ = proc.has_event_; event_.reset( new Event( *proc.event_.get() ) );
      is_point_set_ = false;
      return *this;
    }

    void
    GenericProcess::setPoint( const unsigned int ndim, double* x )
    {
      x_ = std::vector<double>( x, x+ndim );
      is_point_set_ = true;

      if ( CG_LOG_MATCH( "Process:dumpPoint", debugInsideLoop ) )
        dumpPoint();
      clearEvent();
    }

    double
    GenericProcess::x( unsigned int idx ) const
    {
      if ( idx >= x_.size() )
        return -1.;
      return x_[idx];
    }

    void
    GenericProcess::clearEvent()
    {
      event_->restore();
    }

    void
    GenericProcess::setKinematics( const Kinematics& kin )
    {
      kin_ = kin;
      prepareKinematics();
    }

    void
    GenericProcess::prepareKinematics()
    {
      if ( !isKinematicsDefined() )
        throw CG_FATAL( "GenericProcess" ) << "Kinematics not properly defined for the process.";

      const HeavyIon hi1( kin_.incoming_beams.first.pdg ), hi2( kin_.incoming_beams.second.pdg );
      const double m1 = hi1 ? HeavyIon::mass( hi1 ) : PDG::get().mass( kin_.incoming_beams.first.pdg );
      const double m2 = hi2 ? HeavyIon::mass( hi2 ) : PDG::get().mass( kin_.incoming_beams.second.pdg );
      // at some point introduce non head-on colliding beams?
      const auto p1 = Momentum::fromPxPyPzM( 0., 0., +kin_.incoming_beams.first .pz, m1 );
      const auto p2 = Momentum::fromPxPyPzM( 0., 0., -kin_.incoming_beams.second.pz, m2 );
      setIncomingKinematics( p1, p2 );

      s_ = ( p1+p2 ).mass2();
      sqs_ = sqrt( s_ );

      w1_ = p1.mass2();
      w2_ = p2.mass2();

      CG_DEBUG( "GenericProcess" ) << "Kinematics successfully prepared!\n"
        << "  √s = " << sqs_*1.e-3 << " TeV,\n"
        << "  p₁ = " << p1 << ", mass=" << p1.mass() << " GeV\n"
        << "  p₂ = " << p2 << ", mass=" << p2.mass() << " GeV.";
    }

    void
    GenericProcess::dumpPoint() const
    {
      std::ostringstream os;
      for ( unsigned short i = 0; i < x_.size(); ++i )
        os << Form( "  x(%2d) = %8.6f\n\t", i, x_[i] );
      CG_INFO( "GenericProcess" )
        << "Number of integration parameters: " << x_.size() << "\n\t"
        << os.str() << ".";
    }

    void
    GenericProcess::setEventContent( const IncomingState& ini, const OutgoingState& fin )
    {
      if ( !has_event_ )
        return;

      event_->clear();
      //----- add the particles in the event

      //--- incoming state
      for ( const auto& ip : ini ) {
        Particle& p = event_->addParticle( ip.first );
        const auto& part_info = PDG::get()( ip.second );
        p.setPdgId( ip.second, part_info.charge/3. );
        p.setMass( part_info.mass );
        if ( ip.first == Particle::IncomingBeam1
          || ip.first == Particle::IncomingBeam2 )
          p.setStatus( Particle::Status::PrimordialIncoming );
        if ( ip.first == Particle::Parton1
          || ip.first == Particle::Parton2 )
          p.setStatus( Particle::Status::Incoming );
      }
      //--- central system (if not already there)
      const auto& central_system = ini.find( Particle::CentralSystem );
      if ( central_system == ini.end() ) {
        Particle& p = event_->addParticle( Particle::Intermediate );
        p.setPdgId( PDG::invalid );
        p.setStatus( Particle::Status::Propagator );
      }
      //--- outgoing state
      for ( const auto& opl : fin ) { // pair(role, list of PDGids)
        for ( const auto& pdg : opl.second ) {
          Particle& p = event_->addParticle( opl.first );
          const auto& part_info = PDG::get()( pdg );
          p.setPdgId( pdg, part_info.charge/3. );
          p.setMass( part_info.mass );
        }
      }

      //----- define the particles parentage

      const Particles parts = event_->particles();
      for ( const auto& p : parts ) {
        Particle& part = (*event_)[p.id()];
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
      last_event = event_;
    }

    void
    GenericProcess::setIncomingKinematics( const Momentum& p1, const Momentum& p2 )
    {
      if ( !has_event_ || !event_ )
        return;

      CG_DEBUG( "GenericProcess:incomingBeams" )
        << "Incoming primary particles:\n\t"
        << p1 << "\n\t"
        << p2;

      (*event_)[Particle::IncomingBeam1][0].setMomentum( p1 );
      (*event_)[Particle::IncomingBeam2][0].setMomentum( p2 );
    }

    bool
    GenericProcess::isKinematicsDefined()
    {
      if ( !has_event_ )
        return true;

      // check the incoming state
      bool is_incoming_state_set =
        ( !(*event_)[Particle::IncomingBeam1].empty() && !(*event_)[Particle::IncomingBeam2].empty() );

      // check the outgoing state
      bool is_outgoing_state_set =
        ( !(*event_)[Particle::OutgoingBeam1].empty() && !(*event_)[Particle::OutgoingBeam2].empty()
       && !(*event_)[Particle::CentralSystem].empty() );

      // combine both states
      return is_incoming_state_set && is_outgoing_state_set;
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

  //--------------------------------------------------------------------
  // User-friendly display of the kinematics mode
  //--------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const KinematicsMode& pm )
  {
    switch ( pm ) {
      case KinematicsMode::invalid:
        return os << "invalid";
      case KinematicsMode::ElectronElectron:
        return os << "electron/electron";
      case KinematicsMode::ElectronProton:
        return os << "electron/proton";
      case KinematicsMode::ProtonElectron:
        return os << "proton/electron";
      case KinematicsMode::ElasticElastic:
        return os << "elastic/elastic";
      case KinematicsMode::InelasticElastic:
        return os << "inelastic/elastic";
      case KinematicsMode::ElasticInelastic:
        return os << "elastic/inelastic";
      case KinematicsMode::InelasticInelastic:
        return os << "inelastic/inelastic";
    }
    return os;
  }
}
