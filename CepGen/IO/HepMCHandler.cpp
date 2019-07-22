#include "CepGen/IO/ExportHandler.h"
#include "CepGen/IO/HepMCEventInterface.h"

#include "CepGen/Parameters.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Constants.h"

#ifdef HEPMC3
#  include "HepMC3/Version.h"
#  define HEPMC_VERSION HEPMC3_VERSION
#  include "HepMC3/WriterAscii.h"
#  include "HepMC3/WriterAsciiHepMC2.h"
#  include "HepMC3/WriterHEPEVT.h"
#  ifdef HEPMC3_ROOTIO
#    include "HepMC3/WriterRoot.h"
#    include "HepMC3/WriterRootTree.h"
#  endif
using namespace HepMC3;
#else
#  include "HepMC/Version.h"
#  if !defined( HEPMC_VERSION_CODE ) // HepMC v2
#    include "HepMC/IO_GenEvent.h"
#  else
#    include "HepMC/WriterAscii.h"
#    include "HepMC/WriterHEPEVT.h"
#    define HEPMC3
#  endif
using namespace HepMC;
#endif

#include <memory>

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
        void initialise( const Parameters& /*params*/ ) override {}
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
        /// Clear the associated HepMC event content
        void clearEvent();
        /// Populate the associated HepMC event with a Event object
        void fillEvent( const Event& );

        /// Writer object
        std::unique_ptr<T> output_;
        /// Generator cross section and error
        std::shared_ptr<GenCrossSection> xs_;
#ifdef HEPMC3
        /// Auxiliary information on run
        std::shared_ptr<GenRunInfo> runinfo_;
#endif
        /// Associated HepMC event
        std::shared_ptr<CepGenEvent> event_;
    };

    template<typename T>
    HepMCHandler<T>::HepMCHandler( const ParametersList& params ) :
      GenericExportHandler( "hepmc" ),
      output_( new T( params.get<std::string>( "filename", "output.hepmc" ) ) ),
      xs_( new GenCrossSection ),
#ifdef HEPMC3
      runinfo_( new GenRunInfo ),
#endif
      event_( new CepGenEvent )
    {
#ifdef HEPMC3
      output_->set_run_info( runinfo_ );
      runinfo_->set_weight_names( { "Default" } );
#endif
      CG_INFO( "HepMC" )
        << "Interfacing module initialised "
        << "for HepMC version " << HEPMC_VERSION << ".";
    }

    template<typename T> void
    HepMCHandler<T>::operator<<( const Event& evt )
    {
      fillEvent( evt );
#ifdef HEPMC3
      output_->write_event( *event_ );
#else
      output_->write_event( event_.get() );
#endif
    }

    template<typename T> void
    HepMCHandler<T>::setCrossSection( double xsect, double xsect_err )
    {
      xs_->set_cross_section( xsect, xsect_err );
    }

    template<typename T> void
    HepMCHandler<T>::fillEvent( const Event& evt )
    {
      event_->feedEvent( evt );
      // general information
#ifdef HEPMC3
      event_->set_cross_section( xs_ );
      event_->set_run_info( runinfo_ );
#else
      event_->set_cross_section( *xs_ );
#endif
      event_->set_event_number( event_num_++ );
    }
#ifdef HEPMC3
#  if HEPMC_VERSION_CODE >= 3001000
    typedef HepMCHandler<WriterAsciiHepMC2> HepMC2Handler;
#  endif
    typedef HepMCHandler<WriterAscii> HepMC3Handler;
    typedef HepMCHandler<WriterHEPEVT> HEPEVTHandler;
#  ifdef HEPMC3_ROOTIO
    typedef HepMCHandler<WriterRoot> RootHandler;
    typedef HepMCHandler<WriterRootTree> RootTreeHandler;
#  endif
#else
    typedef HepMCHandler<IO_GenEvent> HepMC2Handler;
#endif
  }
}

#ifdef HEPMC3
REGISTER_IO_MODULE( hepmc3, HepMC3Handler )
REGISTER_IO_MODULE( hepmc, HepMC3Handler )
REGISTER_IO_MODULE( hepevt, HEPEVTHandler )
#  ifdef HEPMC3_ROOTIO
REGISTER_IO_MODULE( hepmc_root, RootHandler )
REGISTER_IO_MODULE( hepmc_root_tree, RootTreeHandler )
#  endif
#else
REGISTER_IO_MODULE( hepmc, HepMC2Handler )
#endif
#if defined( HEPMC3 ) && HEPMC_VERSION_CODE >= 3001000
REGISTER_IO_MODULE( hepmc2, HepMC2Handler )
#endif

