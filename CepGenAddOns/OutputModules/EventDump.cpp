#include "CepGen/Modules/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"

#include <iomanip>
#include <fstream>

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Simple event dump module
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jan 2020
     */
    class EventDump : public ExportModule
    {
      public:
        explicit EventDump( const ParametersList& );
        ~EventDump();

        void initialise( const Parameters& ) override;
        void setCrossSection( double, double ) override;
        void operator<<( const Event& ) override;

      private:
        bool save_banner_;
        std::ostream* out_;
    };

    EventDump::EventDump( const ParametersList& params ) :
      ExportModule( params ),
      save_banner_( params.get<bool>( "saveBanner", true ) ),
      out_( nullptr )
    {
      if ( params.has<std::string>( "filename" ) )
        out_ = new std::ofstream( params.get<std::string>( "filename" ) );
      else
        out_ = &std::cout;
    }

    EventDump::~EventDump()
    {
      //file_.close();
    }

    void
    EventDump::initialise( const Parameters& params )
    {
      if ( save_banner_ )
        ( *out_ ) << banner( params, "#" ) << "\n";
    }

    void
    EventDump::setCrossSection( double xsec, double xsec_err )
    {
      ( *out_ ) << "Total cross-section: " << xsec << " +/- " << xsec_err << " pb.\n";
    }

    void
    EventDump::operator<<( const Event& ev )
    {
      ( *out_ ) << ev;
    }
  }
}

REGISTER_IO_MODULE( "dump", EventDump )
