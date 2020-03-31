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
    /// Simple event dump module
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Jan 2020
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
        int print_every_;
        std::ostream* out_;
    };

    EventDump::EventDump( const ParametersList& params ) :
      ExportModule( params ),
      save_banner_( params.get<bool>( "saveBanner", true ) ),
      print_every_( params.get<int>( "printEvery", 1 ) ),
      out_( nullptr )
    {
      if ( params.has<std::string>( "filename" ) )
        out_ = new std::ofstream( params.get<std::string>( "filename" ) );
      else
        out_ = &std::cout;
    }

    EventDump::~EventDump()
    {
      if ( out_ != &std::cout )
        dynamic_cast<std::ofstream*>( out_ )->close();
    }

    void
    EventDump::initialise( const Parameters& params )
    {
      if ( save_banner_ )
        *out_ << banner( params, "#" ) << "\n";
    }

    void
    EventDump::setCrossSection( double xsec, double xsec_err )
    {
      *out_ << "Total cross-section: " << xsec << " +/- " << xsec_err << " pb.\n";
    }

    void
    EventDump::operator<<( const Event& ev )
    {
      if ( print_every_ < 0 || event_num_++ % print_every_ == 0 )
        *out_ << ev << "\n";
    }
  }
}

REGISTER_IO_MODULE( "dump", EventDump )
