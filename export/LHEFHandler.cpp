#include "LHEFHandler.h"

#ifdef HEPMC_VERSION3

namespace OutputHandler
{
  LHEFHandler::LHEFHandler( const char* filename ) :
    HepMCHandler( filename, ExportHandler::LHE )
  {
    output = new LHEF::Writer( filename );
  }

  void
  LHEFHandler::operator<<( const Event* ev )
  {
    fillEvent( ev );
    output->writeEvent();
    clearEvent();
  }
}

#endif
