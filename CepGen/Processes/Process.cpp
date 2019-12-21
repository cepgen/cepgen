#include "CepGen/Processes/Process.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include <iomanip>

namespace cepgen
{
  namespace proc
  {
    Process::Process( const ParametersList& params, const std::string& name, const std::string& description, bool has_event ) :
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      params_( params ), name_( name ), description_( description ),
      first_run( true ), base_jacobian_( 1. ),
      s_( -1. ), sqs_( -1. ),
      mA2_( -1. ), mB2_( -1. ), mX2_( -1. ), mY2_( -1. ),
      t1_( -1. ), t2_( -1. ),
      is_point_set_( false )
    {
      if ( has_event )
        event_.reset( new Event );
    }

    Process::Process( const Process& proc ) :
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      params_( proc.params_ ), name_( proc.name_ ), description_( proc.description_ ),
      first_run( proc.first_run ), base_jacobian_( proc.base_jacobian_ ),
      s_( proc.s_ ), sqs_( proc.sqs_ ),
      mA2_( proc.mA2_ ), mB2_( proc.mB2_ ), mX2_( proc.mX2_ ), mY2_( proc.mY2_ ),
      t1_( -1. ), t2_( -1. ), kin_( proc.kin_ ),
      is_point_set_( false )
    {
      if ( proc.event_ )
        event_.reset( new Event( *proc.event_.get() ) );
    }

    Process&
    Process::operator=( const Process& proc )
    {
      params_ = proc.params_;
      name_ = proc.name_; description_ = proc.description_;
      first_run = proc.first_run;
      base_jacobian_ = proc.base_jacobian_;
      s_ = proc.s_; sqs_ = proc.sqs_;
      mA2_ = proc.mA2_; mB2_ = proc.mB2_;
      mX2_ = proc.mX2_; mY2_ = proc.mY2_;
      kin_ = proc.kin_;
      if ( proc.event_ )
        event_.reset( new Event( *proc.event_.get() ) );
      is_point_set_ = false;
      return *this;
    }

    void
    Process::dumpVariables() const
    {
      std::ostringstream os;
      for ( const auto& var : mapped_variables_ )
        os << "\n\t(" << var.index << ") " << var.type << " mapping (" << var.description << ") in range " << var.limits;
      CG_INFO( "Process:dumpVariables" )
        << "List of variables handled by this kt-factorised process:"
        << os.str();
    }

    Process&
    Process::defineVariable( double& out, const Mapping& type, Limits in, const Limits& default_limits, const std::string& descr )
    {
      if ( !in.valid() ) {
        CG_DEBUG( "Process:defineVariable" )
          << descr << " could not be retrieved from the user configuration!\n\t"
          << "Setting it to the default value: " << default_limits << ".";
        in = default_limits;
      }

      Limits lim = in;
      out = 0.; // reset the variable
      double jacob_weight = 1.; // initialise the local weight for this variable

      switch ( type ) {
        case Mapping::linear:
          jacob_weight = lim.range();
          break;
        case Mapping::square:
          jacob_weight = 2.*lim.range();
          break;
        case Mapping::exponential: {
          lim = { // limits already stored as log(limits)
            ( !lim.hasMin() || lim.min() == 0. ) ? -10. : std::max( log( lim.min() ), -10. ),
            ( !lim.hasMax() || lim.max() == 0. ) ? +10. : std::min( log( lim.max() ), +10. )
          };
          jacob_weight = lim.range(); // use the linear version
        } break;
        case Mapping::power_law:
          jacob_weight = log( lim.max()/lim.min() );
          break;
      }
      mapped_variables_.emplace_back(
        MappingVariable{ descr.empty() ? utils::format( "var%z", mapped_variables_.size() ) : descr,
          lim, out, type, (unsigned short)mapped_variables_.size() } );
      point_coord_.emplace_back( 0. );
      base_jacobian_ *= jacob_weight;
      CG_DEBUG( "Process:defineVariable" )
        << descr << " has been mapped to variable " << mapped_variables_.size() << ".\n\t"
        << "Allowed range for integration: " << in << " (" << lim << ").\n\t"
        << "Variable integration mode: " << type << ".\n\t"
        << "Weight in the Jacobian: " << jacob_weight << ".";
      return *this;
    }

    void
    Process::generateVariables() const
    {
      if ( mapped_variables_.size() == 0 )
        throw CG_FATAL( "Process:vars" )
          << "No variables are mapped for this process!";
      if ( base_jacobian_ == 0. )
        throw CG_FATAL( "Process:vars" )
          << "Point-independant component of the Jacobian for this "
          << "process is null.\n\t"
          << "Please check the validity of the phase space!";

      for ( const auto& var : mapped_variables_ ) {
        if ( !var.limits.valid() )
          continue;
        const double xv = x( var.index ); // between 0 and 1
        switch ( var.type ) {
          case Mapping::linear: {
            var.value = var.limits.x( xv );
          } break;
          case Mapping::exponential: { // limits aleady logarithmic
            var.value = exp( var.limits.x( xv ) ); // transform back to linear
          } break;
          case Mapping::square: {
            var.value = pow( var.limits.x( xv ), 2 );
          } break;
          case Mapping::power_law: {
            const double y = var.limits.max()/var.limits.min();
            var.value = var.limits.min()*pow( y, xv );
          } break;
        }
      }
      if ( CG_LOG_MATCH( "Process:vars", debugInsideLoop ) ) {
        std::ostringstream oss;
        std::string sep;
        for ( const auto& var : mapped_variables_ ) {
          double value = 0.;
          switch ( var.type ) {
            case Mapping::linear:
            case Mapping::exponential:
            case Mapping::power_law:
              value = var.value;
              break;
            case Mapping::square:
              value = sqrt( var.value );
              break;
          }
          oss << sep
            << "variable " << var.index
            << std::left << std::setw( 60 )
            << ( !var.description.empty() ? " ("+var.description+")" : "" )
            << " in range " << std::setw( 20 ) << var.limits
            << " has value " << std::setw( 20 ) << value
            << " (x=" << x( var.index ) << std::right << ")", sep = "\n\t";
        }
        CG_DEBUG_LOOP( "Process:vars" ) << oss.str();
      }
    }

    double
    Process::jacobian() const
    {
      double jac = 1.;
      for ( const auto& var : mapped_variables_ ) {
        if ( !var.limits.valid() )
          continue;
        switch ( var.type ) {
          case Mapping::linear:
            break;
          case Mapping::square:
            jac *= sqrt( var.value );
            break;
          case Mapping::exponential:
          case Mapping::power_law:
            jac *= var.value;
            break;
        }
      }
      return jac;
    }

    void
    Process::setPoint( double* x, const size_t ndim )
    {
      std::copy( x, x+ndim, point_coord_.begin() );
      is_point_set_ = true;

      if ( CG_LOG_MATCH( "Process:dumpPoint", debugInsideLoop ) )
        dumpPoint();
      clearEvent();
    }

    double
    Process::x( unsigned int idx ) const
    {
      try {
        return point_coord_.at( idx );
      } catch ( const std::out_of_range& ) {
        throw CG_FATAL( "Process:x" )
          << "Failed to retrieve coordinate " << idx << " from "
          << "a dimension-" << ndim() << " process!";
      }
    }

    double
    Process::weight()
    {
      if ( !is_point_set_ )
        throw CG_FATAL( "Process:weight" )
          << "Trying to evaluate weight while phase space point\n\t"
          << "coordinates are not set!";

      //--- generate and initialise all variables
      generateVariables();

      //--- compute the integrand
      const double me_integrand = computeWeight();
      //if ( me_integrand <= 0. )
      //  return 0.;

      //--- generate auxiliary (x-dependent) part of the Jacobian for
      //    this phase space point.
      const double aux_jacobian = jacobian();
      if ( aux_jacobian <= 0. )
        return 0.;

      //--- combine every component into a single weight for this point
      const double weight = ( base_jacobian_*aux_jacobian ) * me_integrand;

      CG_DEBUG_LOOP( "Process:weight" )
        << "Jacobian: " << base_jacobian_ << " * " << aux_jacobian
        << " = " << ( base_jacobian_*aux_jacobian ) << ".\n\t"
        << "Integrand = " << me_integrand << "\n\t"
        << "Proc.-specific integrand * Jacobian (excl. global Jacobian) = "
        << ( me_integrand * aux_jacobian ) << "\n\t"
        << "Point weight = " << weight << ".";

      return weight;
    }

    void
    Process::clearEvent()
    {
      event_->restore();
    }

    void
    Process::setKinematics( const Kinematics& kin )
    {
      kin_ = kin;
      //--- initialise the "constant" (wrt x) part of the Jacobian
      base_jacobian_ = 1.;
      mapped_variables_.clear();

      //--- define incoming system
      const HeavyIon hi1( kin_.incoming_beams.first.pdg ), hi2( kin_.incoming_beams.second.pdg );
      const double m1 = hi1 ? HeavyIon::mass( hi1 ) : PDG::get().mass( kin_.incoming_beams.first.pdg );
      const double m2 = hi2 ? HeavyIon::mass( hi2 ) : PDG::get().mass( kin_.incoming_beams.second.pdg );
      // at some point introduce non head-on colliding beams?
      const auto p1 = Momentum::fromPxPyPzM( 0., 0., +kin_.incoming_beams.first .pz, m1 );
      const auto p2 = Momentum::fromPxPyPzM( 0., 0., -kin_.incoming_beams.second.pz, m2 );
      setIncomingKinematics( p1, p2 );

      s_ = ( p1+p2 ).mass2();
      sqs_ = sqrt( s_ );

      mA2_ = p1.mass2();
      mB2_ = p2.mass2();

      CG_DEBUG( "Process" ) << "Kinematics successfully set!\n"
        << "  √s = " << sqs_*1.e-3 << " TeV,\n"
        << "  p1=" << p1 << ",\tmass=" << p1.mass() << " GeV\n"
        << "  p2=" << p2 << ",\tmass=" << p2.mass() << " GeV.";

      //--- process-specific phase space definition
      prepareKinematics();
    }

    void
    Process::dumpPoint() const
    {
      std::ostringstream os;
      for ( unsigned short i = 0; i < point_coord_.size(); ++i )
        os << utils::format( "\n\t  x(%2d) = %8.6f", i, point_coord_[i] );
      CG_INFO( "Process" )
        << "Number of integration parameters: " << mapped_variables_.size()
        << os.str() << ".";
    }

    void
    Process::setEventContent( const IncomingState& ini, const OutgoingState& fin )
    {
      if ( !event_ )
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
        p.setPdgId( (pdgid_t)PDG::invalid );
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
            part.addMother( event_->oneWithRole( Particle::IncomingBeam1 ) );
            break;
          case Particle::OutgoingBeam2:
          case Particle::Parton2:
            part.addMother( event_->oneWithRole( Particle::IncomingBeam2 ) );
            break;
          case Particle::Intermediate:
            part.addMother( event_->oneWithRole( Particle::Parton1 ) );
            part.addMother( event_->oneWithRole( Particle::Parton2 ) );
            break;
          case Particle::CentralSystem:
            part.addMother( event_->oneWithRole( Particle::Intermediate ) );
            break;
          default: break;
        }
      }

      //----- freeze the event as it is

      event_->freeze();
    }

    void
    Process::setIncomingKinematics( const Momentum& p1, const Momentum& p2 )
    {
      if ( !event_ )
        return;

      CG_DEBUG( "Process:incomingBeams" )
        << "Incoming primary particles:\n\t"
        << p1 << "\n\t"
        << p2;

      (*event_)[Particle::IncomingBeam1][0].setMomentum( p1 );
      (*event_)[Particle::IncomingBeam2][0].setMomentum( p2 );
    }

    bool
    Process::isKinematicsDefined()
    {
      if ( !event_ )
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
    operator<<( std::ostream& os, const Process& proc )
    {
      return os << proc.name().c_str();
    }

    std::ostream&
    operator<<( std::ostream& os, const Process* proc )
    {
      return os << proc->name().c_str();
    }

    std::ostream&
    operator<<( std::ostream& os, const Process::Mapping& type )
    {
      switch ( type ) {
        case Process::Mapping::linear: return os << "linear";
        case Process::Mapping::exponential: return os << "exponential";
        case Process::Mapping::square: return os << "squared";
        case Process::Mapping::power_law: return os << "power law";
      }
      return os;
    }
  }
}
