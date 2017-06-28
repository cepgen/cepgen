#ifndef OutputHandler_LHEFHandler_h
#define OutputHandler_LHEFHandler_h

#include "export/HepMCHandler.h"
#include "physics/Event.h"

#ifdef HEPMC_VERSION3

#include "HepMC/LHEF.h"

namespace OutputHandler
{
  /**
   * \brief Handler for the LHE file output
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   * \date Sep 2016
   */
  class LHEFHandler : public HepMCHandler
  {
   public:
    /// Class constructor
    /// \param[in] filename Output file path
    LHEFHandler( const char* filename );
    /// Writer operator
    void operator<<( const Event* );
    
   private:
    /// Writer object (from HepMC)
    LHEF::Writer* output;

  };
}

#endif
#endif
