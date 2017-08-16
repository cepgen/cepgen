#ifndef CepGen_Export_LHEFHandler_h
#define CepGen_Export_LHEFHandler_h

#include "HepMCHandler.h"

#ifdef HEPMC_VERSION3

#include "HepMC/LHEF.h"
#include "CepGen/Physics/Event.h"

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
    void initialise( const Parameters& params );
    /// Writer operator
    void operator<<( const Event* );
    
   private:
    /// Writer object (from HepMC)
    std::unique_ptr<LHEF::Writer> lhe_output_;
    LHEF::HEPRUP run_;
  };
}

#endif
#endif
