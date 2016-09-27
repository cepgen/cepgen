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

    HepMCHandler( const char* );
    ~HepMCHandler();
    void operator<<( const Event* );

   private:

    inline void clearEvent();
    void fillEvent( const Event* );

#ifndef HEPMC_VERSION_CODE
    HepMC::IO_GenEvent* output;
#else
    HepMC::WriterAscii* output;
#endif
    HepMC::GenEvent* event;
    std::vector<HepMC::GenParticle*> particles;
    std::vector<HepMC::GenVertex*> vertices;

  };

}

#endif
