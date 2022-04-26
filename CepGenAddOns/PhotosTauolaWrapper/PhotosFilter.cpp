#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PhotosTauolaWrapper/PhotosTauolaInterface.h"

//#define _LOG_DEBUG_MODE_

#include <Photos/Log.h>
#include <Photos/Photos.h>
#include <Photos/PhotosEvent.h>

using namespace Photospp;

namespace cepgen {
  namespace hadr {
    /// Interface to the Photos decay routine
    class PhotosFilter : public EventModifier {
    public:
      explicit PhotosFilter(const ParametersList&);
      ~PhotosFilter();

      static ParametersDescription description();

      void setRuntimeParameters(const Parameters&) override {}
      void init() override;
      bool run(Event& ev, double& weight, bool full) override;

    private:
      typedef io::PhotosTauolaEvent<PhotosEvent, PhotosParticle> CepGenPhotosEvent;
    };

    PhotosFilter::PhotosFilter(const ParametersList& params) : EventModifier(params) {
      if (steer<bool>("debug"))
        Log::LogAll(true);
      Photos::maxWtInterference(steer<double>("maxWtInterference"));
      Photos::setInfraredCutOff(steer<double>("infraredCutOff"));
      Photos::setInterference(steer<bool>("interference"));
      Photos::setDoubleBrem(steer<bool>("doubleBrem"));
      Photos::setQuatroBrem(steer<bool>("quatroBrem"));
      Photos::setCorrectionWtForW(steer<bool>("correctionWtForW"));
      Photos::setExponentiation(steer<bool>("exponentiation"));
      Photos::setPairEmission(steer<bool>("pairEmission"));
      Photos::setPhotonEmission(steer<bool>("photonEmission"));
      Photos::setMeCorrectionWtForScalar(steer<bool>("meCorrectionWtForScalar"));
      Photos::setMeCorrectionWtForW(steer<bool>("meCorrectionWtForW"));
      Photos::setMeCorrectionWtForZ(steer<bool>("meCorrectionWtForZ"));
      Photos::setTopProcessRadiation(steer<bool>("topProcessRadiation"));
    }

    PhotosFilter::~PhotosFilter() { Log::SummaryAtExit(); }

    void PhotosFilter::init() {
      Photos::setMomentumUnit(Photos::GEV);
      Photos::setAlphaQED(constants::ALPHA_EM);
      //Photos::setSeed( seed_ );
      Photos::initialize();
    }

    bool PhotosFilter::run(Event& ev, double& weight, bool) {
      weight = 1.;

      CepGenPhotosEvent evt(ev, PDG::tau);
      evt.dump();
      evt.process();
      //evt.undecayTaus();
      //evt.decayTaus();
      evt.dump();
      //const auto& pairs = evt[Particle::CentralSystem][0];
      //CG_WARNING("")<<pairs;

      return true;
    }

    ParametersDescription PhotosFilter::description() {
      auto desc = EventModifier::description();
      desc.add<bool>("debug", false).setDescription("log all debugging information?");
      desc.add<double>("maxWtInterference", 1.).setDescription("maximum interference weight");
      desc.add<double>("infraredCutOff", 0.01)
          .setDescription("minimal energy (in units of decaying particle mass) for photons to be explicitly generated");
      desc.add<bool>("interference", true).setDescription("key for interference, matrix element weight");
      desc.add<bool>("doubleBrem", true).setDescription("set double bremsstrahlung generation");
      desc.add<bool>("quatroBrem", false).setDescription("set bremsstrahlung generation up to multiplicity of 4");
      desc.add<bool>("correctionWtForW", true)
          .setDescription("key for partial effects of matrix element (in leptonic W decays)");
      desc.add<bool>("exponentiation", true).setDescription("set exponentiation mode");
      desc.add<bool>("pairEmission", false).setDescription("set pair emission");
      desc.add<bool>("photonEmission", true).setDescription("set photon emission");
      desc.add<bool>("meCorrectionWtForScalar", false)
          .setDescription("switch for complete effects of the matrix element (in scalar to two scalar decays)");
      desc.add<bool>("meCorrectionWtForW", false)
          .setDescription("switch for complete effects of matrix element (in leptonic W decays)");
      desc.add<bool>("meCorrectionWtForZ", false)
          .setDescription("switch for complete effects of matrix element (in leptonic Z decays)");
      desc.add<bool>("topProcessRadiation", true)
          .setDescription("set photon emission in top pair production in quark (gluon) pair annihilation");
      return desc;
    }
  }  // namespace hadr
}  // namespace cepgen

// register event modifier
REGISTER_MODIFIER("photos", PhotosFilter)
