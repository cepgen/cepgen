#ifndef EventWriter_h
#define EventWriter_h

#include "physics/Event.h"
#include "export/ExportHandler.h"

#ifdef HEPMC_LINKED
#include "export/HepMCHandler.h"
#include "export/LHEFHandler.h"
#endif

class EventWriter
{
 public:

  EventWriter( const OutputHandler::ExportHandler::OutputType&, const char* );
  ~EventWriter();

  void SetCrossSection( const float& xsec, const float& err_xsec ) {
#ifdef HEPMC_LINKED
    if ( fFileHandler ) fFileHandler->SetCrossSection( xsec, err_xsec );
#endif
  }
  void SetEventNumber( const unsigned int& ev_id ) {
#ifdef HEPMC_LINKED
    if ( fFileHandler ) fFileHandler->SetEventNumber( ev_id );
#endif
  }

  void operator<<( const Event* );

 private:

  OutputHandler::ExportHandler* fFileHandler;
  OutputHandler::ExportHandler::OutputType fType;

};

#endif
