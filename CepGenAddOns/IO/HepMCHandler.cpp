#include "CepGenAddOns/IO/HepMCEventInterface.h"

#include "CepGen/Parameters.h"

#include "CepGen/Core/ExportHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#ifdef HEPMC3
#  include "HepMC3/Version.h"
#  define HEPMC_VERSION HEPMC3_VERSION
#  define HEPMC_VERSION_CODE HEPMC3_VERSION_CODE
#  include "HepMC3/WriterAscii.h"
#  include "HepMC3/WriterAsciiHepMC2.h"
#  include "HepMC3/WriterHEPEVT.h"
#  ifdef HEPMC3_ROOTIO
#    include "HepMC3/WriterRoot.h"
#    include "HepMC3/WriterRootTree.h"
#  endif
#else
#  include "HepMC/Version.h"
#  if !defined( HEPMC_VERSION_CODE ) // HepMC v2
#    include "HepMC/IO_GenEvent.h"
#    include "HepMC/IO_AsciiParticles.h"
#  else
#    include "HepMC/WriterAscii.h"
#    include "HepMC/WriterHEPEVT.h"
#    define HEPMC3
#  endif
#endif

#include <memory>

using namespace HepMC;

namespace cepgen
{
  namespace io
  {
    /// Handler for the HepMC file output
    /// \tparam T HepMC writer handler (format-dependent)
    /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
    /// \date Sep 2016
    template<typename T>
    class HepMCHandler : public GenericExportHandler
    {
      public:
        /// Class constructor
        explicit HepMCHandler( const ParametersList& );
        ~HepMCHandler();

        void initialise( const Parameters& /*params*/ ) override {}
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
        /// Writer object
        std::unique_ptr<T> output_;
        /// Generator cross section and error
        std::shared_ptr<GenCrossSection> xs_;
#ifdef HEPMC3
        /// Auxiliary information on run
        std::shared_ptr<GenRunInfo> runinfo_;
#endif
    };

    template<typename T>
    HepMCHandler<T>::HepMCHandler( const ParametersList& params ) :
      GenericExportHandler( "hepmc" ),
      output_( new T( params.get<std::string>( "filename", "output.hepmc" ).c_str() ) ),
      xs_( new GenCrossSection )
#ifdef HEPMC3
      , runinfo_( new GenRunInfo )
#endif
    {
#ifdef HEPMC3
      output_->set_run_info( runinfo_ );
      runinfo_->set_weight_names( { "Default" } );
#endif
      CG_INFO( "HepMC" )
        << "Interfacing module initialised "
        << "for HepMC version " << HEPMC_VERSION << ".";
    }

    template<typename T>
    HepMCHandler<T>::~HepMCHandler()
    {
#ifdef HEPMC3
      output_->close();
#endif
    }

    template<typename T> void
    HepMCHandler<T>::operator<<( const Event& evt )
    {
      CepGenEvent event( evt );
      // general information
#ifdef HEPMC3
      event.set_cross_section( xs_ );
      event.set_run_info( runinfo_ );
#else
      event.set_cross_section( *xs_ );
#endif
      event.set_event_number( event_num_++ );
#ifdef HEPMC3
      output_->write_event( event );
#else
      output_->write_event( &event );
#endif
    }

    template<typename T> void
    HepMCHandler<T>::setCrossSection( double xsect, double xsect_err )
    {
      xs_->set_cross_section( xsect, xsect_err );
    }
  }
}

#ifdef HEPMC3
typedef cepgen::io::HepMCHandler<WriterAscii> HepMC3Handler;
typedef cepgen::io::HepMCHandler<WriterHEPEVT> HepMC3HEPEVTHandler;
REGISTER_IO_MODULE( "hepmc", HepMC3Handler )
REGISTER_IO_MODULE( "hepevt", HepMC3HEPEVTHandler )
#  if HEPMC_VERSION_CODE >= 3001000
typedef cepgen::io::HepMCHandler<WriterAsciiHepMC2> HepMC3HepMC2Handler;
REGISTER_IO_MODULE( "hepmc2", HepMC3HepMC2Handler )
#  endif
#  ifdef HEPMC3_ROOTIO
typedef cepgen::io::HepMCHandler<WriterRoot> HepMC3RootHandler;
typedef cepgen::io::HepMCHandler<WriterRootTree> HepMC3RootTreeHandler;
REGISTER_IO_MODULE( "hepmc_root", HepMC3RootHandler )
REGISTER_IO_MODULE( "hepmc_root_tree", HepMC3RootTreeHandler )
#  endif
#else // HepMC version 2 and below
typedef cepgen::io::HepMCHandler<IO_GenEvent> HepMC2Handler;
typedef cepgen::io::HepMCHandler<IO_AsciiParticles> HepMC2AsciiHandler;
REGISTER_IO_MODULE( "hepmc", HepMC2Handler )
REGISTER_IO_MODULE( "hepmc_ascii", HepMC2AsciiHandler )
#endif

