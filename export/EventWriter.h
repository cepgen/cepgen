#ifndef EventWriter_h
#define EventWriter_h

#include "physics/Event.h"

#ifdef HEPMC_LINKED
#include "export/HepMCHandler.h"
#include "export/LHEFHandler.h"
#endif

class EventWriter
{
 public:
  enum OutputType {
    HepMC, LHE
  };

  EventWriter( const OutputType&, const char* );
  ~EventWriter();

  void SetCrossSection( const float& xsec, const float& err_xsec ) {
#ifdef HEPMC_LINKED
    if ( fHepMCHandler ) fHepMCHandler->SetCrossSection( xsec, err_xsec );
    if ( fLHEFHandler ) fLHEFHandler->SetCrossSection( xsec, err_xsec );
#endif
  }
  void SetEventNumber( const unsigned int& ev_id ) {
#ifdef HEPMC_LINKED
    if ( fHepMCHandler ) fHepMCHandler->SetEventNumber( ev_id );
    if ( fLHEFHandler ) fLHEFHandler->SetEventNumber( ev_id );
#endif
  }

  void operator<<( const Event* );

 private:

  OutputType fType;
#ifdef HEPMC_LINKED
  OutputHandler::HepMCHandler* fHepMCHandler;
  OutputHandler::LHEFHandler* fLHEFHandler;
#endif

};

#endif
