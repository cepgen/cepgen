#ifndef Pythia8Hadroniser_h
#define Pythia8Hadroniser_h

#include <algorithm>
#include <Pythia.h>

#include "GenericHadroniser.h"

namespace CepGen
{
  namespace Hadroniser
  {
    /**
     * Full interface to the Pythia8 hadronisation algorithm. It can be used in a single particle decay mode as well as a full event hadronisation using the string model, as in Jetset.
     * \brief Pythia8 hadronisation algorithm
     */
    class Pythia8Hadroniser : public GenericHadroniser
    {
      public:
        Pythia8Hadroniser();
        ~Pythia8Hadroniser();
        bool Hadronise( Particle* part );
        bool Hadronise( Event* ev );

      private:
        bool PrepareHadronisation( Event *ev );
        Pythia8::Pythia* fPy;
    };
  }
}

#endif
