#include "CepGen/Modules/Process.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include <iomanip>

namespace cepgen
{
  namespace proc
  {
    Process::Process( const ParametersList& params, const std::string& name, const std::string& description, bool has_event ) :
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      params_( params ), name_( name ), description_( description ),
      first_run( true ), base_jacobian_( 0. ),
      s_( -1. ), sqs_( -1. ),
      MX_( -1. ), MY_( -1. ), w1_( -1. ), w2_( -1. ),
      t1_( -1. ), t2_( -1. ),
      has_event_( has_event ), event_( new Event ),
      is_point_set_( false )
    {}

    Process::Process( const Process& proc ) :
      mp_( PDG::get().mass( PDG::proton ) ), mp2_( mp_*mp_ ),
      params_( proc.params_ ), name_( proc.name_ ), description_( proc.description_ ),
      first_run( proc.first_run ), base_jacobian_( proc.base_jacobian_ ),
      s_( proc.s_ ), sqs_( proc.sqs_ ),
      MX_( proc.MX_ ), MY_( proc.MY_ ), w1_( proc.w1_ ), w2_( proc.w2_ ),
      t1_( -1. ), t2_( -1. ), kin_( proc.kin_ ),
      has_event_( proc.has_event_ ), event_( new Event( *proc.event_.get() ) ),
      is_point_set_( false )
    {}

    Process&
    Process::operator=( const Process& proc )
    {
      params_ = proc.params_;
      name_ = proc.name_; description_ = proc.description_;
      first_run = proc.first_run;
      base_jacobian_ = proc.base_jacobian_;
      s_ = proc.s_; sqs_ = proc.sqs_;
      MX_ = proc.MX_; MY_ = proc.MY_; w1_ = proc.w1_; w2_ = proc.w2_;
      kin_ = proc.kin_;
      has_event_ = proc.has_event_; event_.reset( new Event( *proc.event_.get() ) );
      is_point_set_ = false;
      return *this;
    }

    Process&
    Process::defineVariable( double& out, const Mapping& type, const Limits& in, Limits default_limits, const char* description )
    {
      Limits lim = in;
      out = 0.; // reset the variable
      if ( !in.valid() ) {
        CG_DEBUG( "Process:defineVariable" )
          << description << " could not be retrieved from the user configuration!\n\t"
          << "Setting it to the default value: " << default_limits << ".";
        lim = default_limits;
      }
      if ( type == Mapping::logarithmic )
        lim = {
          std::max( log( lim.min() ), -10. ),
          std::min( log( lim.max() ), +10. )
        };
      mapped_variables_.emplace_back(
        MappingVariable{ description, lim, out, type, (unsigned short)mapped_variables_.size() } );
      switch ( type ) {
        case Mapping::square:
          base_jacobian_ *= 2.*lim.range();
          break;
        case Mapping::linear: case Mapping::logarithmic:
          base_jacobian_ *= lim.range();
          break;
      }
      CG_DEBUG( "Process:defineVariable" )
        << description << " has been mapped to variable " << mapped_variables_.size() << ".\n\t"
        << "Allowed range for integration: " << lim << ".\n\t"
        << "Variable integration mode: " << type << ".";
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

    double
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

      double jacobian = 1.;
      for ( const auto& cut : mapped_variables_ ) {
        if ( !cut.limits.valid() )
          continue;
        const double xv = x( cut.index ); // between 0 and 1
        switch ( cut.type ) {
          case Mapping::linear: {
            cut.value = cut.limits.x( xv );
          } break;
          case Mapping::logarithmic: {
            cut.value = exp( cut.limits.x( xv ) );
            jacobian *= cut.value;
          } break;
          case Mapping::square: {
            cut.value = cut.limits.x( xv );
            jacobian *= cut.value;
          } break;
        }
      }
      if ( CG_LOG_MATCH( "Process:vars", debugInsideLoop ) ) {
        std::ostringstream oss;
        for ( const auto& cut : mapped_variables_ ) {
          oss << "variable " << cut.index
              << " in range " << std::left << std::setw( 20 ) << cut.limits << std::right
              << " has value " << cut.value << "\n\t";
        }
        CG_DEBUG_LOOP( "Process:vars" ) << oss.str();
      }
      return jacobian;
    }

    void
    Process::setPoint( const unsigned int ndim, double* x )
    {
      for ( size_t i = 0; i < ndim; ++i )
        mapped_variables_[i].value = x[i];
      is_point_set_ = true;

      if ( CG_LOG_MATCH( "Process:dumpPoint", debugInsideLoop ) )
        dumpPoint();
      clearEvent();
    }

    double
    Process::x( unsigned int idx ) const
    {
      if ( idx >= mapped_variables_.size() )
        return -1.;
      return mapped_variables_[idx].value;
    }

    double
    Process::weight()
    {
      //--- generate and initialise all variables, and auxiliary
      //    (x-dependent) part of the Jacobian for this phase space point.
      const double aux_jacobian = generateVariables();
      if ( aux_jacobian <= 0. )
        return 0.;

      //--- compute the integrand and combine together into a single
      //    weight for the phase space point.
      const double me_integrand = computeWeight();
      if ( me_integrand <= 0. )
        return 0.;

      const double weight = ( base_jacobian_*aux_jacobian ) * me_integrand;

      CG_DEBUG_LOOP( "Process:weight" )
        << "Jacobian: " << base_jacobian_ << " * " << aux_jacobian
        << " = " << ( base_jacobian_*aux_jacobian ) << ".\n\t"
        << "Integrand = " << me_integrand << "\n\t"
        << "dW = " << weight << ".";

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

      w1_ = p1.mass2();
      w2_ = p2.mass2();

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
      for ( unsigned short i = 0; i < mapped_variables_.size(); ++i )
        os << Form( "  x(%2d) = %8.6f\n\t", i, mapped_variables_[i].value );
      CG_INFO( "Process" )
        << "Number of integration parameters: " << mapped_variables_.size() << "\n\t"
        << os.str() << ".";
    }

    void
    Process::setEventContent( const IncomingState& ini, const OutgoingState& fin )
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
    Process::setIncomingKinematics( const Momentum& p1, const Momentum& p2 )
    {
      if ( !has_event_ || !event_ )
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
        case Process::Mapping::logarithmic: return os << "logarithmic";
        case Process::Mapping::square: return os << "squared";
      }
      return os;
    }
  }
}
