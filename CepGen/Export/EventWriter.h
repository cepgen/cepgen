#ifndef CepGen_Export_EventWriter_h
#define CepGen_Export_EventWriter_h

#include "CepGen/Parameters.h"
#include "CepGen/Physics/Event.h"
#include "CepGen/Export/ExportHandler.h"

#ifdef HEPMC_LINKED
#include "CepGen/Export/HepMCHandler.h"
#include "CepGen/Export/LHEFHandler.h"
#endif

#include <memory>

namespace CepGen
{
  namespace OutputHandler
  {
    /**
     * \brief Generic events dumper
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class EventWriter
    {
      public:
        /// Class constructor
        /// \param[in] type Requested output type
        /// \param[in] filename Output file path
        EventWriter( const OutputHandler::ExportHandler::OutputType& type, const char* filename );
        ~EventWriter();

        void initialise( const Parameters& params ) { file_handler_->initialise( params ); }

        /// Specify the process cross section and its associated error
        void setCrossSection( const float& xsec, const float& err_xsec ) {
          if ( file_handler_ ) file_handler_->setCrossSection( xsec, err_xsec );
        }
        /// Specify the event number
        void setEventNumber( const unsigned int& ev_id ) {
          if ( file_handler_ ) file_handler_->setEventNumber( ev_id );
        }
        /// Writer operator
        void operator<<( const Event* );

      private:
        /// Inherited file handler
        std::unique_ptr<OutputHandler::ExportHandler> file_handler_;
        /// Type of output requested
        OutputHandler::ExportHandler::OutputType type_;
    };
  }
}

#endif
