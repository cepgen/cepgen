/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <Tauola/Log.h>
#include <Tauola/Tauola.h>
#include <Tauola/TauolaEvent.h>
#include <Tauola/TauolaHepMC3Event.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenHepMC3/CepGenEvent.h"

using namespace Tauolapp;

namespace cepgen::tauola {
  /// Interface to the Tauola decay routine
  class TauolaFilter final : public EventModifier {
  public:
    explicit TauolaFilter(const ParametersList& params)
        : EventModifier(params),
          pol_states_(steer<ParametersList>("polarisations")),
          rad_states_(steer<ParametersList>("radiations")) {
      if (steer<bool>("debug"))
        Log::LogAll(true);
    }
    ~TauolaFilter() override { Log::SummaryAtExit(); }

    void initialise() override {
      Tauola::setUnits(Tauola::GEV, Tauola::MM);
      Tauola::initialize();
      Tauola::setSeed(seed_, 2. * seed_, 4. * seed_);
      Tauola::momentum_conservation_threshold = 1.e-6;
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
      const auto rad_cutoff = rad_states_.get<double>("cutoff", 0.01);
      if (rad_cutoff > 0.)
        Tauola::setRadiationCutOff(rad_cutoff);  // default energy is 0.01 (in units of half the decaying particle mass)

      //--- default parameters
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
    bool run(Event& event, double& weight, bool /*fast*/) override {
      weight = 1.;
      HepMC3::CepGenEvent hepmc_event(event);
      TauolaHepMC3Event tauola_event(&hepmc_event);
      tauola_event.decayTaus();
      hepmc_event.dump();
      hepmc_event.merge(event);
      return true;
    }

    static ParametersDescription description() {
      auto desc = EventModifier::description();
      desc.setDescription("Tauola interface");
      desc.add("debug", false).setDescription("debugging mode");

      auto pol_desc = ParametersDescription();
      pol_desc.add("full", true);
      pol_desc.add("GAMMA", Tauola::spin_correlation.GAMMA);
      pol_desc.add("Z0", Tauola::spin_correlation.Z0);
      pol_desc.add("HIGGS", Tauola::spin_correlation.HIGGS);
      pol_desc.add("HIGGS_H", Tauola::spin_correlation.HIGGS_H);
      pol_desc.add("HIGGS_A", Tauola::spin_correlation.HIGGS_A);
      pol_desc.add("HIGGS_PLUS", Tauola::spin_correlation.HIGGS_PLUS);
      pol_desc.add("HIGGS_MINUS", Tauola::spin_correlation.HIGGS_MINUS);
      pol_desc.add("W_PLUS", Tauola::spin_correlation.W_PLUS);
      pol_desc.add("W_MINUS", Tauola::spin_correlation.W_MINUS);
      desc.add("polarisations", pol_desc);

      auto rad_desc = ParametersDescription();
      rad_desc.add("enable", false);
      rad_desc.add("cutoff", -1.);
      desc.add("radiations", rad_desc);

      desc.add("sameParticleDecayMode", -1);
      desc.add("oppositeParticleDecayMode", -1);

      auto br_desc = ParametersDescription();
      br_desc.add("mode", -1).setDescription("decay mode");
      br_desc.add("branchingRatio", 0.).setDescription("branching fraction");
      desc.addParametersDescriptionVector("branchingRatios", br_desc, {});
      return desc;
    }

  private:
    const ParametersList pol_states_, rad_states_;
  };
}  // namespace cepgen::tauola
using cepgen::tauola::TauolaFilter;
REGISTER_MODIFIER("tauola", TauolaFilter);  // register event modifier
