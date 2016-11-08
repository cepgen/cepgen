#include "LHEFHandler.h"

#if HEPMC_VERSION_CODE>=3000000

namespace OutputHandler
{
  LHEFHandler::LHEFHandler( const char* filename ) :
    ExportHandler( ExportHandler::LHE )
  {
    output = new LHEF::Writer(filename);
  }

  void
  LHEFHandler::operator<<( const Event* ev )
  {
    fillEvent(ev);
    output->writeEvent();
    clearEvent();
  }

  void
  LHEFHandler::fillEvent( const Event* ev )
  {
    // ... do whatever is needed for output->hepeup
    ConstParticlesRef part_vec = ev->GetConstParticlesRef();
    //HEPEUT*
    for ( unsigned int i=0; i<part_vec.size(); i++ ) {

    }
  }

  void
  LHEFHandler::clearEvent()
  {
    // ...
    output->hepeup.clear();
  }
}

#endif
