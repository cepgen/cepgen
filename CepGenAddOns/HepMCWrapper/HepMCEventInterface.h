#ifndef CepGenAddOns_EventInterfaces_HepMCEventInterface_h
#define CepGenAddOns_EventInterfaces_HepMCEventInterface_h

#ifdef HEPMC3
#include "HepMC3/GenEvent.h"
#define HepMC HepMC3
#else
#include "HepMC/GenEvent.h"
#endif
#include <unordered_map>

namespace cepgen {
  class Event;
}

namespace HepMC {
  /// Interfacing between CepGen and HepMC event definitions
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class CepGenEvent : public GenEvent {
  public:
    /// Construct an event interface from a CepGen Event object
    CepGenEvent(const cepgen::Event& ev);

  private:
    std::unordered_map<unsigned short, std::shared_ptr<GenParticle> > assoc_map_;
  };
}  // namespace HepMC
#endif
