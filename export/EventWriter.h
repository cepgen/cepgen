#ifndef EventWriter_h
#define EventWriter_h

#include "physics/Event.h"

#include "export/HepMCHandler.h"
#include "export/LHEFHandler.h"

class EventWriter
{
 public:
  enum OutputType {
    HepMC, LHE
  };

  EventWriter( const OutputType&, const char* );
  ~EventWriter();

  void SetCrossSection( const float& xsec, const float& err_xsec ) {
    if ( fHepMCHandler ) fHepMCHandler->SetCrossSection( xsec, err_xsec );
  }
  void SetEventNumber( const unsigned int& ev_id ) {
    if ( fHepMCHandler ) fHepMCHandler->SetEventNumber( ev_id );
  }

  void operator<<( const Event* );

 private:

  OutputType fType;
  OutputHandler::HepMCHandler* fHepMCHandler;
  OutputHandler::LHEFHandler* fLHEFHandler;

};

#endif
