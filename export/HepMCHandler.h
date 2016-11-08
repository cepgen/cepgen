#ifndef OutputHandler_HepMCHandler_h
#define OutputHandler_HepMCHandler_h

#include "export/ExportHandler.h"

#include "HepMC/Version.h"

#ifndef HEPMC_VERSION_CODE
#include "HepMC/IO_GenEvent.h"
#else
#include "HepMC/WriterAscii.h"
#endif

#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/GenParticle.h"

namespace OutputHandler
{
  /**
   * \brief Handler for the HepMC file output
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   * \date Sep 2016
   */
  class HepMCHandler : public ExportHandler
  {
   public:
    /// Class constructor
    /// \param[in] filename Output file path
    HepMCHandler( const char* filename );
    ~HepMCHandler();
    /// Writer operator
    void operator<<( const Event* );

   private:
    /// Clear the associated HepMC event content
    inline void clearEvent();
    /// Populate the associated HepMC event with a Event object
    void fillEvent( const Event* );

#ifndef HEPMC_VERSION_CODE
    /// Writer object (from HepMC v>=3)
    HepMC::IO_GenEvent* output;
#else
    /// Writer object (from HepMC v<3)
    HepMC::WriterAscii* output;
#endif
    /// Associated HepMC event
    HepMC::GenEvent* event;
    /// List of particles in the event
    std::vector<HepMC::GenParticle*> particles;
    /// List of vertices in the event
    std::vector<HepMC::GenVertex*> vertices;

  };

}

#endif
