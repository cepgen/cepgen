#ifndef EventWriter_h
#define EventWriter_h

#include "physics/Event.h"
#include "export/ExportHandler.h"

#ifdef HEPMC_LINKED
#include "export/HepMCHandler.h"
#include "export/LHEFHandler.h"
#endif

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

    /// Specify the process cross section and its associated error
    void setCrossSection( const float& xsec, const float& err_xsec ) {
#ifdef HEPMC_LINKED
      if ( file_handler_ ) file_handler_->setCrossSection( xsec, err_xsec );
#endif
    }
    /// Specify the event number
    void setEventNumber( const unsigned int& ev_id ) {
#ifdef HEPMC_LINKED
      if ( file_handler_ ) file_handler_->setEventNumber( ev_id );
#endif
    }
    /// Writer operator
    void operator<<( const Event* );

   private:

    /// Inherited file handler
    OutputHandler::ExportHandler* file_handler_;
    /// Type of output requested
    OutputHandler::ExportHandler::OutputType type_;

  };
}

#endif
