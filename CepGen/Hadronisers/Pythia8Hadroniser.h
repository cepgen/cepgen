#ifndef Pythia8Hadroniser_h
#define Pythia8Hadroniser_h

#include "GenericHadroniser.h"

#ifdef PYTHIA8
#include <Pythia8/Pythia.h>
#include <memory>
#endif

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
        bool hadronise( const Particle& part, Event& ev ) override;
        bool hadronise( Event& ev ) override;

        void setSeed( long long seed );
        void readString( const char* param ) {
#ifdef PYTHIA8
          pythia_->readString( param );
#endif
        }
        void readString( const std::string& param ) {
#ifdef PYTHIA8
          pythia_->readString( param );
#endif
        }
        bool init() { return pythia_->init(); }

      private:
#ifdef PYTHIA8
        std::unique_ptr<Pythia8::Pythia> pythia_;
#endif
    };
  }
}

#endif
