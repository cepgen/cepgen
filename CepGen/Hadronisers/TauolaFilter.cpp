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
        const ParametersList pol_states_, rad_states_;
        std::unique_ptr<HepMC::CepGenEvent> event_;
    };

    TauolaFilter::TauolaFilter( const ParametersList& params ) :
      GenericHadroniser( params, "tauola" ),
      pol_states_( params.get<ParametersList>( "polarisations" ) ),
      rad_states_( params.get<ParametersList>( "radiations" ) ),
      event_( new HepMC::CepGenEvent )
    {}

    void
    TauolaFilter::init()
    {
      Tauola::setUnits( Tauola::GEV, Tauola::CM );
      //Tauola::setSeed( seed_ );
      Tauola::initialize();
      //--- spin correlations
      Tauola::spin_correlation.setAll( pol_states_.get<bool>( "all", true ) );
      Tauola::spin_correlation.GAMMA = pol_states_.get<bool>( "gamma", true );
      Tauola::spin_correlation.Z0 = pol_states_.get<bool>( "Z0", true );
      Tauola::spin_correlation.HIGGS = pol_states_.get<bool>( "H", true );
      Tauola::spin_correlation.HIGGS_H = pol_states_.get<bool>( "H_H", true );
      Tauola::spin_correlation.HIGGS_A = pol_states_.get<bool>( "H_A", true );
      Tauola::spin_correlation.HIGGS_PLUS = pol_states_.get<bool>( "H_p", true );
      Tauola::spin_correlation.HIGGS_MINUS = pol_states_.get<bool>( "H_m", true );
      Tauola::spin_correlation.W_PLUS = pol_states_.get<bool>( "W_p", true );
      Tauola::spin_correlation.W_MINUS = pol_states_.get<bool>( "W_m", true );
      //--- radiation states
      Tauola::setRadiation( rad_states_.get<bool>( "enable", true ) );
      const auto rad_cutoff = rad_states_.get<double>( "cutoff", -1. );
      if ( rad_cutoff > 0. )
        Tauola::setRadiationCutOff( rad_cutoff );
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

