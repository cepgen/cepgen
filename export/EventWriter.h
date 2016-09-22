#ifndef EventWriter_h
#define EventWriter_h

#include "../include/Event.h"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/GenParticle.h"

/*#include "LHE/..."*/

class EventWriter
{
 public:
  enum OutputType {
    HepMC=0/*, LHE*/
  };

  EventWriter(const OutputType&, const char*);
  ~EventWriter();

  void SetCrossSection(const float& xsec, const float& err_xsec) {
    fCrossSect = xsec;
    fCrossSectErr = err_xsec;
  }
  void SetEventNumber(const unsigned int& ev_id) { fEventNum = ev_id; }

  void operator<<(const Event*);

 private:
  HepMC::GenParticle getHepMCParticle(const Particle*) const;
  HepMC::GenEvent getHepMCEvent(const Event*) const;

  OutputType fType;

  HepMC::IO_GenEvent* fHepMCOutput;

  float fCrossSect, fCrossSectErr;
  unsigned int fEventNum;
};

#endif
