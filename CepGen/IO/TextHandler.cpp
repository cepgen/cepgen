#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include <fstream>
#include <regex>

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class TextHandler : public GenericExportHandler
    {
      public:
        explicit TextHandler( const ParametersList& );
        ~TextHandler();

        void initialise( const Parameters& ) override;
        void operator<<( const Event& ) override;

      private:
        /// Retrieve a named variable from a particle
        double variable( const Particle&, const std::string& ) const;
        /// Retrieve a named variable from the whole event
        double variable( const Event&, const std::string& ) const;

        static const std::regex rgx_select_id_, rgx_select_role_;
        static constexpr double INVALID_OUTPUT = -999.;

        std::ofstream file_;
        const std::vector<std::string> variables_;
        const bool print_banner_, print_variables_;
        const std::string separator_;

        //--- variables definition
        typedef std::pair<unsigned short,std::string> IndexedVariable;
        std::unordered_map<short,std::vector<IndexedVariable> > variables_per_id_;
        std::unordered_map<Particle::Role,std::vector<IndexedVariable> > variables_per_role_;
        std::vector<IndexedVariable> variables_for_event_;
        unsigned short num_vars_;
        std::ostringstream oss_vars_;

        //--- auxiliary helper maps
        const std::unordered_map<std::string,Particle::Role> role_str_ = {
          { "ib1", Particle::Role::IncomingBeam1 }, { "ib2", Particle::Role::IncomingBeam2 },
          { "ob1", Particle::Role::OutgoingBeam1 }, { "ob2", Particle::Role::OutgoingBeam2 },
          { "pa1", Particle::Role::Parton1 }, { "pa2", Particle::Role::Parton2 },
          { "cs",  Particle::Role::CentralSystem },
          { "int", Particle::Role::Intermediate }
        };
        typedef double( Particle::Momentum::*pMethod )(void) const;
        /// Mapping of string variables to momentum getter methods
        const std::unordered_map<std::string,pMethod> m_mom_str_ = {
          { "px",  &Particle::Momentum::px },
          { "py",  &Particle::Momentum::py },
          { "pz",  &Particle::Momentum::pz },
          { "pt",  &Particle::Momentum::pt },
          { "eta", &Particle::Momentum::eta },
          { "phi", &Particle::Momentum::phi },
          { "m",   &Particle::Momentum::mass },
          { "e",   &Particle::Momentum::energy },
          { "p",   &Particle::Momentum::p },
          { "pt2", &Particle::Momentum::pt2 },
          { "th",  &Particle::Momentum::theta },
          { "y",   &Particle::Momentum::rapidity }
        };

        //--- kinematic variables
        double sqrts_;
        unsigned long num_evts_;
    };

    const std::regex TextHandler::rgx_select_id_( "(\\w+)\\((\\d+)\\)" );
    const std::regex TextHandler::rgx_select_role_( "(\\w+)\\(([a-z]+\\d?)\\)" );

    TextHandler::TextHandler( const ParametersList& params ) :
      GenericExportHandler( "text" ),
      file_           ( params.get<std::string>( "filename", "output.txt" ) ),
      variables_      ( params.get<std::vector<std::string> >( "variables" ) ),
      print_banner_   ( params.get<bool>( "saveBanner", true ) ),
      print_variables_( params.get<bool>( "saveVariables", true ) ),
      separator_      ( params.get<std::string>( "separator", "\t" ) ),
      num_vars_( 0 )
    {
      std::smatch sm;
      oss_vars_.clear();
      std::string sep;
      for ( const auto& var : variables_ ) {
        if ( std::regex_match( var, sm, rgx_select_id_ ) )
          variables_per_id_[stod( sm[2].str() )].emplace_back( std::make_pair( num_vars_, sm[1].str() ) );
        else if ( std::regex_match( var, sm, rgx_select_role_ ) ) {
          const auto& str_role = sm[2].str();
          if ( role_str_.count( str_role ) == 0 ) {
            CG_WARNING( "TextHandler" )
              << "Invalid particle role retrieved from configuration: \"" << str_role << "\".\n\t"
              << "Skipping the variable \"" << var << "\" in the output module.";
            continue;
          }
          variables_per_role_[role_str_.at( str_role )].emplace_back( std::make_pair( num_vars_, sm[1].str() ) );
        }
        else // event-level variables
          variables_for_event_.emplace_back( std::make_pair( num_vars_, var ) );
        oss_vars_ << sep << var, sep = separator_;
        ++num_vars_;
      }
    }

    TextHandler::~TextHandler()
    {
      file_.close();
    }

    void
    TextHandler::initialise( const Parameters& params )
    {
      sqrts_ = params.kinematics.sqrtS();
      num_evts_ = 0ul;
      if ( print_banner_ )
        file_ << banner( params, "#" ) << "\n";
      if ( print_variables_ )
        file_ << "# " << oss_vars_.str() << "\n";
    }

    void
    TextHandler::operator<<( const Event& ev )
    {
      std::vector<double> vars( num_vars_ );
      //--- extract and order the variables to be retrieved
      //--- particle-level variables (indexed by integer id)
      for ( const auto& id_vars : variables_per_id_ ) {
        const auto& part = ev[id_vars.first];
        //--- loop over the list of variables for this particle
        for ( const auto& var : id_vars.second )
          vars[var.first] = variable( part, var.second );
      }
      //--- particle-level variables (indexed by role)
      for ( const auto& role_vars : variables_per_role_ ) {
        const auto& part = ev[role_vars.first][0];
        //--- loop over the list of variables for this particle
        for ( const auto& var : role_vars.second )
          vars[var.first] = variable( part, var.second );
      }
      //--- event-level variables
      for ( const auto& var : variables_for_event_ )
        vars[var.first] = variable( ev, var.second );
      //--- write down the variables list in the file
      std::string sep;
      for ( const auto& var : vars )
        file_ << sep << var, sep = separator_;
      file_ << "\n";
      ++num_evts_;
    }

    double
    TextHandler::variable( const Particle& part, const std::string& var ) const
    {
      if ( m_mom_str_.count( var ) ) {
        auto meth = m_mom_str_.at( var );
        return ( part.momentum().*meth )();
      }
      if ( var == "xi"  ) return 1.-part.energy()*2./sqrts_;
      if ( var == "pdg" ) return (double)part.integerPdgId();
      if ( var == "charge" ) return part.charge();
      if ( var == "status" ) return (double)part.status();
      CG_WARNING( "TextHandler" )
        << "Failed to retrieve variable \"" << var << "\".";
      return INVALID_OUTPUT;
    }

    double
    TextHandler::variable( const Event& ev, const std::string& var ) const
    {
      if ( var == "np" )
        return (double)ev.size();
      if ( var == "nev" )
        return (double)num_evts_+1;
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
      CG_WARNING( "TextHandler" )
        << "Failed to retrieve the event-level variable \"" << var << "\".";
      return INVALID_OUTPUT;
    }
  }
}

REGISTER_IO_MODULE( text, TextHandler )
