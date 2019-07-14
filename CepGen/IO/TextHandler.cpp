#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"

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

        static const std::regex rgx_select_id_, rgx_select_role_;
        static constexpr double INVALID_OUTPUT = -999.;

        std::ofstream file_;
        std::vector<std::string> variables_;
        std::unordered_map<short,std::vector<std::pair<unsigned short,std::string> > > variables_per_id_;
        std::unordered_map<Particle::Role,std::vector<std::pair<unsigned short,std::string> > > variables_per_role_;
        unsigned short num_vars_;

        //--- kinematic variables
        double sqrts_;
    };

    const std::regex TextHandler::rgx_select_id_( "(\\w+)\\((\\d+)\\)" );
    const std::regex TextHandler::rgx_select_role_( "(\\w+)\\(([a-z]+\\d?)\\)" );

    TextHandler::TextHandler( const ParametersList& params ) :
      GenericExportHandler( "text" ),
      file_( params.get<std::string>( "filename", "output.txt" ) ),
      variables_( params.get<std::vector<std::string> >( "variables" ) ),
      num_vars_( 0 )
    {
      const auto vars_tmp = variables_;
      variables_.clear();
      std::smatch sm;
      std::string sep;
      file_ << "# ";
      for ( const auto& var : vars_tmp ) {
        if ( std::regex_match( var, sm, rgx_select_id_ ) )
          variables_per_id_[stod( sm[2].str() )].emplace_back( std::make_pair( num_vars_, sm[1].str() ) );
        else if ( std::regex_match( var, sm, rgx_select_role_ ) ) {
          const auto& str_role = sm[2].str();
          auto role = Particle::Role::UnknownRole;
          if      ( str_role == "ib1" ) role = Particle::Role::IncomingBeam1;
          else if ( str_role == "ib2" ) role = Particle::Role::IncomingBeam2;
          else if ( str_role == "ob1" ) role = Particle::Role::OutgoingBeam1;
          else if ( str_role == "ob2" ) role = Particle::Role::OutgoingBeam2;
          else if ( str_role == "cs"  ) role = Particle::Role::CentralSystem;
          else if ( str_role == "int" ) role = Particle::Role::Intermediate;
          else if ( str_role == "pa1" ) role = Particle::Role::Parton1;
          else if ( str_role == "pa2" ) role = Particle::Role::Parton2;
          else {
            CG_WARNING( "TextHandler" )
              << "Invalid particle role retrieved from configuration: \"" << str_role << "\".\n\t"
              << "Skipping the variable \"" << var << "\" in the output module.";
            continue;
          }
          variables_per_role_[role].emplace_back( std::make_pair( num_vars_, sm[1].str() ) );
        }
        else {
          CG_WARNING( "TextHandler" )
            << "Generic variables retrieval not yet supported.\n\t"
            << "Skipping the variable \"" << var << "\" in the output module.";
          variables_.emplace_back( var );
        }
        file_ << sep << var, sep = "\t";
        ++num_vars_;
      }
      file_ << "\n";
    }

    TextHandler::~TextHandler()
    {
      file_.close();
    }

    void
    TextHandler::initialise( const Parameters& params )
    {
      sqrts_ = params.kinematics.sqrtS();
    }

    void
    TextHandler::operator<<( const Event& ev )
    {
      std::vector<double> vars( num_vars_ );
      //--- extract and order the variables to be retrieved
      for ( const auto& id_vars : variables_per_id_ ) {
        //--- first get the particle
        const auto& part = ev[id_vars.first];
        //--- then loop on the variables
        for ( const auto& var : id_vars.second )
          vars[var.first] = variable( part, var.second );
      }
      for ( const auto& role_vars : variables_per_role_ ) {
        //--- first get the particle
        const auto& part = ev[role_vars.first][0];
        //--- then loop on the variables
        for ( const auto& var : role_vars.second )
          vars[var.first] = variable( part, var.second );
      }
      //--- write down the variables list in the file
      std::string separator;
      for ( const auto& var : vars )
        file_ << separator << var, separator = "\t";
      file_ << "\n";
    }

    double
    TextHandler::variable( const Particle& part, const std::string& var ) const
    {
      if      ( var == "px"  ) return part.momentum().px();
      else if ( var == "py"  ) return part.momentum().py();
      else if ( var == "pz"  ) return part.momentum().pz();
      else if ( var == "pt"  ) return part.momentum().pt();
      else if ( var == "m"   ) return part.mass();
      else if ( var == "e"   ) return part.energy();
      else if ( var == "xi"  ) return 1.-part.energy()*2./sqrts_;
      else if ( var == "eta" ) return part.momentum().eta();
      else if ( var == "phi" ) return part.momentum().phi();
      else if ( var == "status" ) return (double)part.status();
      CG_WARNING( "TextHandler" )
        << "Failed to retrieve variable \"" << var << "\".";
      return INVALID_OUTPUT;
    }
  }
}

REGISTER_IO_MODULE( text, TextHandler )
