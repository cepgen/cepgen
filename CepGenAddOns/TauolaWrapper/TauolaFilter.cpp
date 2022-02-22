#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/TauolaWrapper/PhotosTauolaInterface.h"

//#define _LOG_DEBUG_MODE_

#include <Tauola/Log.h>
#include <Tauola/Tauola.h>
#include <Tauola/TauolaEvent.h>

using namespace Tauolapp;

namespace cepgen {
  namespace hadr {
    /// Interface to the Tauola decay routine
    class TauolaFilter : public EventModifier {
    public:
      explicit TauolaFilter(const ParametersList&);
      ~TauolaFilter();

      void setRuntimeParameters(const Parameters&) override {}
      void init() override;
      bool run(Event& ev, double& weight, bool full) override;

      void setCrossSection(double /* xsec */, double /* xsec_err */) override {}

      static ParametersDescription description();

    private:
      const ParametersList pol_states_, rad_states_;
      typedef io::PhotosTauolaEvent<TauolaEvent, TauolaParticle> CepGenTauolaEvent;
    };

    TauolaFilter::TauolaFilter(const ParametersList& params)
        : EventModifier(params),
          pol_states_(steer<ParametersList>("polarisations")),
          rad_states_(steer<ParametersList>("radiations")) {
      Log::LogAll(true);
    }

    TauolaFilter::~TauolaFilter() { Log::SummaryAtExit(); }

    void TauolaFilter::init() {
      Tauola::setUnits(Tauola::GEV, Tauola::MM);
      Tauola::initialize();
      Tauola::setSeed(seed_, 2. * seed_, 4. * seed_);
      if (!Tauola::getIsTauolaIni())
        throw CG_FATAL("TauolaFilter:init") << "Tauola was not properly initialised!";

      //--- spin correlations
      if (pol_states_.has<bool>("full"))
        Tauola::spin_correlation.setAll(pol_states_.get<bool>("full"));
      pol_states_.fill<bool>("GAMMA", Tauola::spin_correlation.GAMMA);
      pol_states_.fill<bool>("Z0", Tauola::spin_correlation.Z0);
      pol_states_.fill<bool>("HIGGS", Tauola::spin_correlation.HIGGS);
      pol_states_.fill<bool>("HIGGS_H", Tauola::spin_correlation.HIGGS_H);
      pol_states_.fill<bool>("HIGGS_A", Tauola::spin_correlation.HIGGS_A);
      pol_states_.fill<bool>("HIGGS_PLUS", Tauola::spin_correlation.HIGGS_PLUS);
      pol_states_.fill<bool>("HIGGS_MINUS", Tauola::spin_correlation.HIGGS_MINUS);
      pol_states_.fill<bool>("W_PLUS", Tauola::spin_correlation.W_PLUS);
      pol_states_.fill<bool>("W_MINUS", Tauola::spin_correlation.W_MINUS);

      //--- radiation states
      if (rad_states_.has<bool>("enable"))
        Tauola::setRadiation(rad_states_.get<bool>("enable"));
      const auto rad_cutoff = rad_states_.get<double>("cutoff", -1.);
      if (rad_cutoff > 0.)
        Tauola::setRadiationCutOff(rad_cutoff);  // default energy is 0.01 (in units of half the decaying particle mass)

      //--- default parameters
      //Tauola::setDecayingParticle( 15 );
      Tauola::setSameParticleDecayMode(steer<int>("sameParticleDecayMode"));
      Tauola::setOppositeParticleDecayMode(steer<int>("oppositeParticleDecayMode"));

      //--- list of tau decay branching fractions
      for (const auto& br_per_mode : steer<std::vector<ParametersList> >("branchingRatios")) {
        const auto mode = br_per_mode.get<int>("mode");
        const auto br = br_per_mode.get<double>("branchingRatio");
        Tauola::setTauBr(mode, br);
        CG_DEBUG("TauolaFilter:init") << "Branching ratio for mode " << mode << " set to " << br << ".";
      }
    }

    bool TauolaFilter::run(Event& ev, double& weight, bool /* full */) {
      weight = 1.;

      CepGenTauolaEvent evt(ev, PDG::tau);
      evt.dump();
      //evt.undecayTaus();
      //evt.decayTaus();
      //evt.dump();
      for (auto* tau : evt.findParticles(PDG::tau))
        Tauola::decayOne(tau);
      throw CG_FATAL("") << "fini";
      //const auto& pairs = evt[Particle::CentralSystem][0];
      //CG_WARNING("")<<pairs;

      return true;
    }

    ParametersDescription TauolaFilter::description() {
      auto desc = EventModifier::description();
      desc.setDescription("Tauola interface");

      auto pol_desc = ParametersDescription();
      pol_desc.add<bool>("full", true);
      pol_desc.add<bool>("GAMMA", Tauola::spin_correlation.GAMMA);
      pol_desc.add<bool>("Z0", Tauola::spin_correlation.Z0);
      pol_desc.add<bool>("HIGGS", Tauola::spin_correlation.HIGGS);
      pol_desc.add<bool>("HIGGS_H", Tauola::spin_correlation.HIGGS_H);
      pol_desc.add<bool>("HIGGS_A", Tauola::spin_correlation.HIGGS_A);
      pol_desc.add<bool>("HIGGS_PLUS", Tauola::spin_correlation.HIGGS_PLUS);
      pol_desc.add<bool>("HIGGS_MINUS", Tauola::spin_correlation.HIGGS_MINUS);
      pol_desc.add<bool>("W_PLUS", Tauola::spin_correlation.W_PLUS);
      pol_desc.add<bool>("W_MINUS", Tauola::spin_correlation.W_MINUS);
      desc.add<ParametersDescription>("polarisations", pol_desc);

      auto rad_desc = ParametersDescription();
      rad_desc.add<bool>("enable", false);
      rad_desc.add<double>("cutoff", -1.);
      desc.add<ParametersDescription>("radiations", rad_desc);

      desc.add<int>("sameParticleDecayMode", -1);
      desc.add<int>("oppositeParticleDecayMode", -1);

      auto br_desc = ParametersDescription();
      br_desc.add<int>("mode", -1).setDescription("decay mode");
      br_desc.add<double>("branchingRatio", 0.).setDescription("branching fraction");
      desc.addParametersDescriptionVector("branchingRatios", br_desc, {});
      return desc;
    }
  }  // namespace hadr
}  // namespace cepgen

// register event modifier
REGISTER_MODIFIER("tauola", TauolaFilter)
