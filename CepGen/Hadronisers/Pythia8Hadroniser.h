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

        bool hadronise( Event& ev, double& weight ) override;
        void setSeed( long long seed );

#ifdef PYTHIA8
        double decay( const Particle& part, Event& ev );
        bool fragment( const Particle& part, Event& ev );

        void readString( const char* param ) {
          pythia_->readString( param );
        }
        void readString( const std::string& param ) {
          pythia_->readString( param );
        }
        bool init() {
          return pythia_->init();
        }
#endif

      private:
#ifdef PYTHIA8
        std::unique_ptr<Pythia8::Pythia> pythia_;
#endif
    };
  }
}

#endif
