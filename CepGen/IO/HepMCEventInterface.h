#ifndef CepGen_IO_HepMCEventInterface_h
#define CepGen_IO_HepMCEventInterface_h

#ifdef HEPMC3
# include "HepMC3/GenEvent.h"
# define HepMC HepMC3
#else
# include "HepMC/GenEvent.h"
#endif
#include <unordered_map>

namespace cepgen { class Event; }

namespace HepMC
{
  /// Interfacing between CepGen and HepMC event definitions
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class CepGenEvent : public GenEvent
  {
    public:
      CepGenEvent( const cepgen::Event& ev );

    private:
#ifdef HEPMC3
      std::unordered_map<unsigned short,std::shared_ptr<GenParticle> > assoc_map_;
#else
      std::unordered_map<unsigned short,GenParticle*> assoc_map_;
#endif
  };
}
#endif
