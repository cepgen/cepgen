#ifndef CepGen_IO_HepMCEventInterface_h
#define CepGen_IO_HepMCEventInterface_h

#ifdef HEPMC3
#  include "HepMC3/GenEvent.h"
#  define HepMC HepMC3
#else
#  include "HepMC/GenEvent.h"
#endif


namespace cepgen { class Event; }

namespace HepMC
{
  /// Interfacing between CepGen and HepMC event definitions
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class CepGenEvent : public GenEvent
  {
    public:
      explicit CepGenEvent();
      /// Feed a new CepGen event to this conversion object
      /// \param[in] ev CepGen event to be fed
      void feedEvent( const cepgen::Event& ev );
  };
}
#endif

