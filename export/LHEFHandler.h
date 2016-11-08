#ifndef OutputHandler_LHEFHandler_h
#define OutputHandler_LHEFHandler_h

#include "export/ExportHandler.h"
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
  class LHEFHandler : public ExportHandler
  {
   public:
    /// Class constructor
    /// \param[in] filename Output file path
    LHEFHandler( const char* filename );
    /// Writer operator
    void operator<<( const Event* );
    
   private:
    /// Fill the handler with the original Event object
    void fillEvent( const Event* );
    /// Remove all references to the original Event object
    void clearEvent();
    /// Writer object (from HepMC)
    LHEF::Writer* output;

  };
}

#endif
#endif
