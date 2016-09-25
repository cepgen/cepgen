#ifndef OutputHandler_LHEFHandler_h
#define OutputHandler_LHEFHandler_h

#ifdef HEPMC_VERSION_CODE && HEPMC_VERSION_CODE>=3000000

#include "HepMC/LHEF.h"

namespace OutputHandler
{
  class LHEFHandler
  {
   public:

    LHEFHandler(const char* filename) : HepMCHandler(filename) {;}
    
   private:

    //LHEF

  };
}

#endif
#endif
