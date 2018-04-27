#include "CepGen/Processes/GenericProcess.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Physics/PDG.h"

namespace CepGen
{
  namespace Process
  {
    const double GenericProcess::mp_ = ParticleProperties::mass( PDG::Proton );
    const double GenericProcess::mp2_ = GenericProcess::mp_*GenericProcess::mp_;

    GenericProcess::GenericProcess( const std::string& name, const std::string& description, bool has_event ) :
      name_( name ), description_( description ),
      first_run( true ),
      s_( -1. ), sqs_( -1. ),
      MX_( -1. ), MY_( -1. ), w1_( -1. ), w2_( -1. ),
      t1_( -1. ), t2_( -1. ),
      has_event_( has_event ), event_( new Event ),
      is_point_set_( false )
    {}

    GenericProcess::GenericProcess( const GenericProcess& proc ) :
      name_( proc.name_ ), description_( proc.description_ ),
      first_run( proc.first_run ),
      s_( proc.s_ ), sqs_( proc.sqs_ ),
      MX_( proc.MX_ ), MY_( proc.MY_ ), w1_( proc.w1_ ), w2_( proc.w2_ ),
      t1_( -1. ), t2_( -1. ), cuts_( proc.cuts_ ),
      has_event_( proc.has_event_ ), event_( new Event( *proc.event_.get() ) ),
      is_point_set_( false )
    {}

    GenericProcess::~GenericProcess()
    {}

    GenericProcess&
    GenericProcess::operator=( const GenericProcess& proc )
    {
      name_ = proc.name_; description_ = proc.description_;
      first_run = proc.first_run;
      s_ = proc.s_; sqs_ = proc.sqs_;
      MX_ = proc.MX_; MY_ = proc.MY_; w1_ = proc.w1_; w2_ = proc.w2_;
      cuts_ = proc.cuts_;
      has_event_ = proc.has_event_; event_.reset( new Event( *proc.event_.get() ) );
      is_point_set_ = false;
      return *this;
    }

    void
    GenericProcess::setPoint( const unsigned int ndim, double* x )
    {
      x_ = std::vector<double>( x, x+ndim );
      is_point_set_ = true;

      if ( CG_EXCEPT_MATCH( "Process:dumpPoint", debugInsideLoop ) )
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
        throw CG_FATAL( "GenericProcess" ) << "Kinematics not properly defined for the process.";

      // at some point introduce non head-on colliding beams?
      Particle::Momentum p1( 0., 0.,  cuts_.inp.first ), p2( 0., 0., -cuts_.inp.second );
      // on-shell beam particles
      p1.setMass( ParticleProperties::mass( cuts_.inpdg.first ) );
      p2.setMass( ParticleProperties::mass( cuts_.inpdg.second ) );
      setIncomingKinematics( p1, p2 );

      sqs_ = CMEnergy( p1, p2 );
      s_ = sqs_*sqs_;

      w1_ = p1.mass2();
      w2_ = p2.mass2();

      CG_DEBUG( "GenericProcess" ) << "Kinematics successfully prepared! sqrt(s) = " << sqs_ << ".";
    }

    void
    GenericProcess::dumpPoint() const
    {
      std::ostringstream os;
      for ( unsigned short i = 0; i < x_.size(); ++i ) {
        os << Form( "  x(%2d) = %8.6f\n\t", i, x_[i] );
      }
      CG_INFO( "GenericProcess" )
        << "Number of integration parameters: " << x_.size() << "\n\t"
        << os.str() << ".";
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
        p.setPdgId( PDG::invalid );
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
        case Kinematics::Mode::ElectronElectron: {
          fp1 = FormFactors::Trivial(); // electron (trivial) form factor
          fp2 = FormFactors::Trivial(); // electron (trivial) form factor
        } break;
        case Kinematics::Mode::ProtonElectron: {
          fp1 = FormFactors::ProtonElastic( -t1_ ); // proton elastic form factor
          fp2 = FormFactors::Trivial(); // electron (trivial) form factor
        } break;
        case Kinematics::Mode::ElectronProton: {
          fp1 = FormFactors::Trivial(); // electron (trivial) form factor
          fp2 = FormFactors::ProtonElastic( -t2_ ); // proton elastic form factor
        } break;
        case Kinematics::Mode::ElasticElastic: {
          fp1 = FormFactors::ProtonElastic( -t1_ ); // proton elastic form factor
          fp2 = FormFactors::ProtonElastic( -t2_ ); // proton elastic form factor
        } break;
        case Kinematics::Mode::ElasticInelastic: {
          fp1 = FormFactors::ProtonElastic( -t1_ );
          fp2 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t2_, w2_, my2 );
        } break;
        case Kinematics::Mode::InelasticElastic: {
          fp1 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t1_, w1_, mx2 );
          fp2 = FormFactors::ProtonElastic( -t2_ );
        } break;
        case Kinematics::Mode::InelasticInelastic: {
          fp1 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t1_, w1_, mx2 );
          fp2 = FormFactors::ProtonInelastic( cuts_.structure_functions, -t2_, w2_, my2 );
        } break;
      }
    }

    bool
    GenericProcess::isKinematicsDefined()
    {
      // check the incoming state
      bool is_incoming_state_set =
        ( !event_->getByRole( Particle::IncomingBeam1 ).empty()
       && !event_->getByRole( Particle::IncomingBeam2 ).empty() );

      // check the outgoing state
      bool is_outgoing_state_set =
        ( !event_->getByRole( Particle::OutgoingBeam1 ).empty()
       && !event_->getByRole( Particle::OutgoingBeam2 ).empty()
       && !event_->getByRole( Particle::CentralSystem ).empty() );

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
}
