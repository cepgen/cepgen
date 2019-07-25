#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

namespace cepgen
{
  namespace utils
  {
    const std::regex EventBrowser::rgx_select_id_( "(\\w+)\\((\\d+)\\)" );
    const std::regex EventBrowser::rgx_select_role_( "(\\w+)\\(([a-z]+\\d?)\\)" );

    double
    EventBrowser::get( const Event& ev, const std::string& var ) const
    {
      std::smatch sm;
      if ( std::regex_match( var, sm, rgx_select_id_ ) ) { // per-id variable
        const auto& var_name = sm[1].str();
        const auto& part = ev[std::stod( sm[2].str() )];
        return variable( part, var_name );
      }
      else if ( std::regex_match( var, sm, rgx_select_role_ ) ) { // per-role variable
        const auto& var_name = sm[1].str();
        const auto& str_role = sm[2].str();
        if ( role_str_.count( str_role ) == 0 ) {
          CG_WARNING( "TextHandler" )
            << "Invalid particle role retrieved from configuration: \"" << str_role << "\".\n\t"
            << "Skipping the variable \"" << var << "\" in the output module.";
          return INVALID_OUTPUT;
        }
        const auto& part = ev[role_str_.at( str_role )][0];
        return variable( part, var_name );
      }
      else // event-level variable
        return variable( ev, var );
    }

    double
    EventBrowser::variable( const Particle& part, const std::string& var ) const
    {
      if ( m_mom_str_.count( var ) ) {
        auto meth = m_mom_str_.at( var );
        return ( part.momentum().*meth )();
      }
      //if ( var == "xi"  ) return 1.-part.momentum().energy()*2./sqrts_;
      if ( var == "pdg" ) return (double)part.integerPdgId();
      if ( var == "charge" ) return part.charge();
      if ( var == "status" ) return (double)part.status();
      CG_WARNING( "EventBrowser" )
        << "Failed to retrieve variable \"" << var << "\".";
      return INVALID_OUTPUT;
    }

    double
    EventBrowser::variable( const Event& ev, const std::string& var ) const
    {
      if ( var == "np" )
        return (double)ev.size();
      //if ( var == "nev" )
      //  return (double)num_evts_+1;
      if ( var == "nob1" || var == "nob2" ) {
        unsigned short out = 0.;
        for ( const auto& part : ev[
          var == "nob1"
          ? Particle::Role::OutgoingBeam1
          : Particle::Role::OutgoingBeam2
        ] )
          if ( (int)part.status() > 0 )
            out++;
        return (double)out;
      }
      if ( var == "tgen" )
        return ev.time_generation;
      if ( var == "ttot" )
        return ev.time_total;
      CG_WARNING( "EventBrowser" )
        << "Failed to retrieve the event-level variable \"" << var << "\".";
      return INVALID_OUTPUT;
    }
  }
}

