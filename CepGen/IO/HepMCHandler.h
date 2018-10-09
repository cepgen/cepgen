#ifndef CepGen_IO_HepMCHandler_h
#define CepGen_IO_HepMCHandler_h

#include "CepGen/IO/ExportHandler.h"

#ifdef LIBHEPMC
#  include "HepMC/Version.h"
#  ifndef HEPMC_VERSION_CODE // HepMC v2
#    include "HepMC/IO_GenEvent.h"
#    include "HepMC/SimpleVector.h"
#  else // HepMC v3+
#    define HEPMC_VERSION3
#    include "HepMC/WriterAscii.h"
#    include "HepMC/FourVector.h"
#  endif
#  include "HepMC/GenEvent.h"
#endif

#include <memory>

namespace CepGen
{
  namespace output
  {
    /**
     * \brief Handler for the HepMC file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class HepMCHandler : public ExportHandler
    {
      public:
        /// Class constructor
        /// \param[in] filename Output file path
        /// \param[in] type Output type
        HepMCHandler( const char* filename, const ExportHandler::OutputType& type = ExportHandler::HepMC );
        void initialise( const Parameters& params ) override {}
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
        /// Clear the associated HepMC event content
        void clearEvent();
        /// Populate the associated HepMC event with a Event object
        void fillEvent( const Event& );
#ifdef LIBHEPMC
#  ifdef HEPMC_VERSION3
        /// Writer object (from HepMC v3+)
        std::unique_ptr<HepMC::WriterAscii> output_;
        /// Generator cross section and error
        HepMC::GenCrossSectionPtr xs_;
#  else
        /// Writer object (from HepMC v<3)
        std::unique_ptr<HepMC::IO_GenEvent> output_;
        /// Generator cross section and error
        HepMC::GenCrossSection xs_;
#  endif
        /// Associated HepMC event
        std::shared_ptr<HepMC::GenEvent> event_;
#endif
    };
  }
}

#endif
