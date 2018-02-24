#ifndef CepGen_Export_LHEFHandler_h
#define CepGen_Export_LHEFHandler_h

#ifndef LIBHEPMC
#error "HepMC is not linked to this instance!"
#endif

#include "HepMCHandler.h"

#ifndef HEPMC_VERSION3
#error "HepMC v3 is required for the LHEF export!"
#else

#include "ExportHandler.h"
#include "HepMC/LHEF.h"
#include "CepGen/Event/Event.h"

namespace CepGen
{
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
      explicit LHEFHandler( const char* filename );
      void initialise( const Parameters& params );
      /// Writer operator
      void operator<<( const Event* );

     private:
      /// Writer object (from HepMC)
      std::unique_ptr<LHEF::Writer> lhe_output_;
      LHEF::HEPRUP run_;
    };
  }
}

#endif
#endif
