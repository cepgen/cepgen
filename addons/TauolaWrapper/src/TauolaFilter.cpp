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
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/StreamCollector.h"
#include "CepGen/Utils/String.h"
#include "CepGenHepMC3/CepGenEvent.h"

using namespace Tauolapp;

namespace cepgen::tauola {
  std::unique_ptr<utils::RandomGenerator> gRandomGenerator;  ///< CepGen-specific random number generator to use

  /// Interface to the Tauola decay routine
  class TauolaFilter final : public EventModifier {
  public:
    explicit TauolaFilter(const ParametersList& params) : EventModifier(params) {
      if (const auto& random_generator = steer<ParametersList>("randomGenerator"); !random_generator.empty()) {
        gRandomGenerator = RandomGeneratorFactory::get().build(random_generator);
        Tauola::setRandomGenerator([]() -> double { return gRandomGenerator->uniform(); });
      }
      Log::LogAll(steer<bool>("debug"));
      std::string buf;
      {
        Tauola::setUnits(Tauola::GEV, Tauola::MM);
        auto sc = utils::StreamCollector(buf);
        Tauola::initialize();
      }
      CG_INFO("TauolaFilter") << "Tauola initialised. Output:\n" << std::string(80, '-') << buf << std::string(80, '-');
      if (!Tauola::getIsTauolaIni())
        throw CG_FATAL("TauolaFilter") << "Tauola was not properly initialised!";
      // default parameters
      Tauola::setSeed(seed_, 2. * seed_, 4. * seed_);
      Tauola::momentum_conservation_threshold = steer<double>("momentumConservationThreshold");
      Tauola::setDecayingParticle(steer<int>("decayingParticle"));
      Tauola::setSameParticleDecayMode(steer<int>("sameParticleDecayMode"));
      Tauola::setOppositeParticleDecayMode(steer<int>("oppositeParticleDecayMode"));
      // list of polarisation and spin correlations-specific parameters
      if (const auto& pol_states = steer<ParametersList>("polarisations"); !pol_states.empty()) {  // spin correlations
        if (pol_states.has<bool>("full"))
          Tauola::spin_correlation.setAll(pol_states.get<bool>("full"));
        pol_states.fill("GAMMA", Tauola::spin_correlation.GAMMA);
        pol_states.fill("Z0", Tauola::spin_correlation.Z0);
        pol_states.fill("HIGGS", Tauola::spin_correlation.HIGGS);
        pol_states.fill("HIGGS_H", Tauola::spin_correlation.HIGGS_H);
        pol_states.fill("HIGGS_A", Tauola::spin_correlation.HIGGS_A);
        pol_states.fill("HIGGS_PLUS", Tauola::spin_correlation.HIGGS_PLUS);
        pol_states.fill("HIGGS_MINUS", Tauola::spin_correlation.HIGGS_MINUS);
        pol_states.fill("W_PLUS", Tauola::spin_correlation.W_PLUS);
        pol_states.fill("W_MINUS", Tauola::spin_correlation.W_MINUS);
      }
      // list of enabled radiation states
      if (const auto& rad_states = steer<ParametersList>("radiations"); !rad_states.empty()) {
        if (rad_states.has<bool>("enable"))
          Tauola::setRadiation(rad_states.get<bool>("enable"));
        if (const auto rad_cutoff = rad_states.get<double>("cutoff", 0.01); rad_cutoff > 0.)
          Tauola::setRadiationCutOff(rad_cutoff);
      }
      // list of tau decay branching fractions
      for (const auto& br_per_mode : steer<std::vector<ParametersList> >("branchingRatios")) {
        const auto mode = br_per_mode.get<int>("mode");
        if (const auto br = br_per_mode.get<double>("branchingRatio"); br > 0.) {
          Tauola::setTauBr(mode, br);
          CG_DEBUG("TauolaFilter") << "Branching ratio for mode " << mode << " set to " << br << ".";
        }
      }
    }
    ~TauolaFilter() override { Log::SummaryAtExit(); }

    static ParametersDescription description() {
      auto desc = EventModifier::description();
      desc.setDescription("Tauola interface");
      desc.add("debug", false).setDescription("debugging mode");
      desc.add("decayingParticle", 15).setDescription("pdg id of the particle to decay (+-15 typically)");
      desc.add("sameParticleDecayMode", static_cast<int>(Tauola::All))
          .allow(Tauola::All, "all")
          .allow(Tauola::ElectronMode, "electron")
          .allow(Tauola::MuonMode, "muon")
          .allow(Tauola::PionMode, "pion")
          .allow(Tauola::RhoMode, "rho")
          .allow(Tauola::A1Mode, "A_1")
          .allow(Tauola::KMode, "K")
          .allow(Tauola::KStarMode, "K*")
          .setDescription("uniformise the decay mode of all particle with the one given in 'decayingParticle'");
      desc.add("oppositeParticleDecayMode", static_cast<int>(Tauola::All))
          .allow(Tauola::All, "all")
          .allow(Tauola::ElectronMode, "electron")
          .allow(Tauola::MuonMode, "muon")
          .allow(Tauola::PionMode, "pion")
          .allow(Tauola::RhoMode, "rho")
          .allow(Tauola::A1Mode, "A_1")
          .allow(Tauola::KMode, "K")
          .allow(Tauola::KStarMode, "K*")
          .setDescription(
              "uniformise the decay mode of all particle with opposite charge to the one given in 'decayingParticle'");
      desc.add("momentumConservationThreshold", 1.e-6)
          .setDescription("numerical limit to ensure momentum conservation");
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
      rad_desc.add("enable", false).setDescription("switch on/off bremsstrahlung in leptonic tau decays?");
      rad_desc.add("cutoff", -1.)
          .setDescription(
              "radiation energy cut-off above which photon is explicitly generated (in units of half the decaying "
              "particle mass)");
      desc.add("radiations", rad_desc).setDescription("Bremsstrahlung parameters block");
      auto br_desc = ParametersDescription();
      br_desc.add("mode", -1).setDescription("decay mode");
      br_desc.add("branchingRatio", 0.).setDescription("branching fraction");
      desc.addParametersDescriptionVector("branchingRatios", br_desc, {})
          .setDescription("List of decay-specific branching fractions");
      desc.add("randomGenerator", ParametersDescription{}).setDescription("overridden random generator algorithm");
      return desc;
    }

    bool run(Event& event, double& weight, bool /*fast*/) override {
      weight = 1.;
      CepGenEvent hepmc_event(event);  // conversion to a HepMC3 format
      const auto hepmc_event_size_before = hepmc_event.particles_size();
      TauolaHepMC3Event tauola_event(&hepmc_event);
      tauola_event.decayTaus();
      hepmc_event.merge(event);  // merge everything back into the original event
      if (hepmc_event.particles_size() == hepmc_event_size_before)
        return false;
      return true;
    }
  };
}  // namespace cepgen::tauola
using cepgen::tauola::TauolaFilter;
REGISTER_MODIFIER("tauola", TauolaFilter);  // register event modifier
