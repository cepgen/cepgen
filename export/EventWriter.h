#ifndef EventWriter_h
#define EventWriter_h

#include "../include/Event.h"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/GenParticle.h"

namespace EventWriter
{
  namespace HepMC
  {
    typedef ::HepMC::IO_GenEvent output;
    ::HepMC::GenEvent* Event(const Event&);

    static ::HepMC::GenEvent* last_event = 0;
  }
}

/*#include "LHE/..."

namespace EventWriter
{
  namespace LHE
  {
  }
  }*/

#endif
