#ifndef OutputHandler_LHEFHandler_h
#define OutputHandler_LHEFHandler_h

#include "export/HepMCHandler.h"
#include "physics/Event.h"

#ifdef HEPMC_VERSION3

#include "HepMC/LHEF.h"
#include "HepMC/HEPEVT_Wrapper.h"

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
    HEPEVT hepevt_buf_;
    /// Writer object (from HepMC)
    std::unique_ptr<LHEF::Writer> lhe_output_;

  };
}

#endif
#endif
