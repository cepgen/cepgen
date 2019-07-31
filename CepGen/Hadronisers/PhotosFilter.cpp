#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/EventModifierHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/IO/PhotosTauolaInterface.h"

//#define _LOG_DEBUG_MODE_

#include <Photos/PhotosEvent.h>
#include <Photos/Photos.h>
#include <Photos/Log.h>

using namespace Photospp;

namespace cepgen
{
  namespace hadr
  {
    /// Interface to the Photos decay routine
    class PhotosFilter : public EventModifier
    {
      public:
        explicit PhotosFilter( const ParametersList& );
        ~PhotosFilter();

        void setParameters( const Parameters& ) override {}
        void init() override;
        bool run( Event& ev, double& weight, bool full ) override;

        void setCrossSection( double xsec, double xsec_err ) override {}

      private:
        typedef io::PhotosTauolaEvent<PhotosEvent,PhotosParticle> CepGenPhotosEvent;
    };

    PhotosFilter::PhotosFilter( const ParametersList& params ) :
      EventModifier( params, "photos" )
    {
      Log::LogAll( true );
      if ( params.has<double>( "maxWtInterference" ) )
        // maximum interference weight
        Photos::maxWtInterference( params.get<double>( "maxWtInterference" ) );
      if ( params.has<double>( "infraredCutOff" ) )
        // minimal energy (in units of decaying particle mass) for photons to be explicitly generated
        Photos::setInfraredCutOff( params.get<double>( "infraredCutOff" ) );
      if ( params.has<bool>( "interference" ) )
        // key for interference, matrix element weight
        Photos::setInterference( params.get<bool>( "interference" ) );
      if ( params.has<bool>( "doubleBrem" ) )
        // set double bremsstrahlung generation
        Photos::setDoubleBrem( params.get<bool>( "doubleBrem" ) );
      if ( params.has<bool>( "quatroBrem" ) )
        // set bremsstrahlung generation up to multiplicity of 4
        Photos::setQuatroBrem( params.get<bool>( "quatroBrem" ) );
      if ( params.has<bool>( "correctionWtForW" ) )
        // key for partial effects of  matrix element (in leptonic W decays)
        Photos::setCorrectionWtForW( params.get<bool>( "correctionWtForW" ) );
      if ( params.has<bool>( "exponentiation" ) )
        // set exponentiation mode
        Photos::setExponentiation( params.get<bool>( "exponentiation" ) );
      if ( params.has<bool>( "pairEmission" ) )
        // set pair emission
        Photos::setPairEmission( params.get<bool>( "pairEmission" ) );
      if ( params.has<bool>( "photonEmission" ) )
        // set photon emission
        Photos::setPhotonEmission( params.get<bool>( "photonEmission" ) );
      if ( params.has<bool>( "meCorrectionWtForScalar" ) )
        // switch for complete effects of matrix element (in  scalar  to 2 scalars decays)
        Photos::setMeCorrectionWtForScalar( params.get<bool>( "meCorrectionWtForScalar" ) );
      if ( params.has<bool>( "meCorrectionWtForW" ) )
        // switch for complete effects of matrix element (in leptonic W decays)
        Photos::setMeCorrectionWtForW( params.get<bool>( "meCorrectionWtForW" ) );
      if ( params.has<bool>( "meCorrectionWtForZ" ) )
        // switch for complete effects of matrix element (in leptonic Z decays)
        Photos::setMeCorrectionWtForZ( params.get<bool>( "meCorrectionWtForZ" ) );
      if ( params.has<bool>( "topProcessRadiation" ) )
        // set photon emission in top pair production in quark (gluon) pair annihilation
        Photos::setTopProcessRadiation( params.get<bool>( "topProcessRadiation" ) );
    }

    PhotosFilter::~PhotosFilter()
    {
      Log::SummaryAtExit();
    }

    void
    PhotosFilter::init()
    {
      Photos::setMomentumUnit( Photos::GEV );
      Photos::setAlphaQED( constants::ALPHA_EM );
      //Photos::setSeed( seed_ );
      Photos::initialize();
    }

    bool
    PhotosFilter::run( Event& ev, double& weight, bool full )
    {
      weight = 1.;

      CepGenPhotosEvent evt( ev, PDG::tau );
      evt.dump();
      //evt.undecayTaus();
      //evt.decayTaus();
      evt.dump();
      //const auto& pairs = evt[Particle::CentralSystem][0];
      //CG_WARNING("")<<pairs;

      return true;
    }
  }
}

// register event modifier
REGISTER_MODIFIER( "photos", PhotosFilter )

