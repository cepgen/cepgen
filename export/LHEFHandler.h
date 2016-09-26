#ifndef OutputHandler_LHEFHandler_h
#define OutputHandler_LHEFHandler_h

#include "physics/Event.h"

#include "HepMC/Version.h"

#if HEPMC_VERSION_CODE>=3000000

#include "HepMC/LHEF.h"

namespace OutputHandler
{
  /**
   * \brief Handler for the LHE file output
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   * \date Sep 2016
   */
  class LHEFHandler
  {
   public:

    LHEFHandler( const char* filename );
    void operator<<( const Event* );
    
   private:

    void fillEvent( const Event* );
    void clearEvent();

    LHEF::Writer* output;

  };
}

#endif
#endif
