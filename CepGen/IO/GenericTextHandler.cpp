#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

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
    class GenericTextHandler : public GenericExportHandler
    {
      public:
        explicit GenericTextHandler( const ParametersList& );
        ~GenericTextHandler();

        void initialise( const Parameters& ) override;
        void operator<<( const Event& ) override;

      private:
        static const std::regex rgx_select_id_, rgx_select_role_;
        std::ofstream file_;
        std::vector<std::string> variables_;
        std::unordered_map<short,std::string> variables_per_role_, variables_per_id_;
    };

    const std::regex GenericTextHandler::rgx_select_id_( "\\w\\((\\d+)\\)" );
    const std::regex GenericTextHandler::rgx_select_role_( "\\w\\(([a-z]+\\d?)\\)" );

    GenericTextHandler::GenericTextHandler( const ParametersList& params ) :
      GenericExportHandler( "text" ),
      file_( params.get<std::string>( "filename", "output.txt" ) ),
      variables_( params.get<std::vector<std::string> >( "variables" ) )
    {
      const auto vars_tmp = variables_;
      variables_.clear();
      std::smatch sm;
      for ( const auto& var : vars_tmp ) {
        if ( std::regex_match( var, sm, rgx_select_id_ ) )
          CG_INFO("")<<sm.str();
      }
    }

    GenericTextHandler::~GenericTextHandler()
    {
      file_.close();
    }

    void
    GenericTextHandler::initialise( const Parameters& params )
    {
    }

    void
    GenericTextHandler::operator<<( const Event& ev )
    {
    }
  }
}

REGISTER_IO_MODULE( text, GenericTextHandler )
