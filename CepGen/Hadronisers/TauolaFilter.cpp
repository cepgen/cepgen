#include "CepGen/Hadronisers/GenericHadroniser.h"
#include "CepGen/Hadronisers/HadronisersHandler.h"

#include "CepGen/IO/HepMCEventInterface.h"

#include "CepGen/Core/ParametersList.h"

#include <Tauola/Tauola.h>
#include <Tauola/TauolaHepMCEvent.h>

using namespace Tauolapp;

namespace cepgen
{
  namespace hadr
  {
    /// Interface to the Tauola decay routine
    class TauolaFilter : public GenericHadroniser
    {
      public:
        explicit TauolaFilter( const ParametersList& );

        void setParameters( const Parameters& ) override {}
        inline void readString( const char* param ) override {}
        void init() override;
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}

      private:
        std::unique_ptr<HepMC::CepGenEvent> event_;
    };

    TauolaFilter::TauolaFilter( const ParametersList& plist ) :
      GenericHadroniser( plist, "tauola" ),
      event_( new HepMC::CepGenEvent )
    {}

    void
    TauolaFilter::init()
    {
      Tauola::setUnits( Tauola::GEV, Tauola::CM );
      //Tauola::setSeed( seed_ );
      Tauola::initialize();
    }

    bool
    TauolaFilter::run( Event& ev, double& weight, bool full )
    {
      weight = 1.;

      event_->feedEvent( ev );
      TauolaHepMCEvent evt( event_.get() );
      //evt.undecayTaus();
      evt.decayTaus();

      return true;
    }
  }
}

// register hadroniser and define alias
REGISTER_HADRONISER( tauola, TauolaFilter )

