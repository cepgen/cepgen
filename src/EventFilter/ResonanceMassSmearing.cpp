/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;

/// À-posteriori Breit-Wigner smearing of the masses of an on-shell resonance pair
class ResonanceMassSmearing final : public EventModifier {
public:
  explicit ResonanceMassSmearing(const ParametersList& params)
      : EventModifier(params),
        pdg_ids_(steer<std::vector<int>>("pdgIds")),
        num_widths_(steer<double>("numWidths")),
        random_generator_(RandomGeneratorFactory::get().build(steer<ParametersList>("randomGenerator"))) {}

  static ParametersDescription description() {
    auto desc = EventModifier::description();
    desc.setDescription("Breit-Wigner mass smearing of a pair of on-shell central resonances");
    desc.add("pdgIds", std::vector<int>{24}).setDescription("absolute PDG ids of the resonances to be smeared");
    desc.add("numWidths", 10.).setDescription("truncation of the mass range, in units of the resonance width");
    desc.add("maxTrials", 100).setDescription("maximal number of mass sampling attempts per event");
    desc.add("randomGenerator", RandomGeneratorFactory::get().describeParameters("stl"))
        .setDescription("random number generator engine");
    return desc;
  }

  bool run(Event& event, double& weight, bool /*fast*/) override {
    weight = 1.;
    if (!event.hasRole(Particle::Role::CentralSystem))
      return true;
    std::vector<Particle*> pair;
    for (auto& part_ref : event[Particle::Role::CentralSystem])
      if (auto& part = part_ref.get();
          std::find(pdg_ids_.begin(), pdg_ids_.end(), std::abs(static_cast<int>(part.integerPdgId()))) !=
          pdg_ids_.end())
        pair.emplace_back(&part);
    if (pair.size() != 2)  // only (pair-wise) resonances are handled
      return true;
    if (!pair.at(0)->children().empty() || !pair.at(1)->children().empty()) {
      CG_WARNING("ResonanceMassSmearing:run") << "Resonances already decayed: mass smearing must precede the decay.";
      return true;  // would otherwise break the parentage momentum balance
    }

    const auto pdg_id = pair.at(0)->pdgId();
    const auto m0 = PDG::get().mass(pdg_id), width = PDG::get().width(pdg_id);
    if (width <= 0.)
      return true;

    Momentum p_tot;  // pair four-momentum and invariant mass
    for (const auto* part : pair)
      p_tot += part->momentum();
    const auto mass_pair2 = p_tot.mass2();
    if (mass_pair2 <= 0.)
      return false;
    const auto mass_pair = std::sqrt(mass_pair2);

    // truncated Breit-Wigner masses, with the pair kept above threshold
    const auto m_min = std::max(0., m0 - num_widths_ * width), m_max = m0 + num_widths_ * width;
    double m1{0.}, m2{0.};
    bool found = false;
    for (int trial = 0; trial < max_trials_ && !found; ++trial) {
      m1 = random_generator_->breitWigner(m0, 0.5 * width);
      m2 = random_generator_->breitWigner(m0, 0.5 * width);
      found = m1 > m_min && m1 < m_max && m2 > m_min && m2 < m_max && m1 + m2 < mass_pair;
    }
    if (!found)
      return false;

    // lab -> pair rest frame
    const auto p_tot_inv = Momentum::fromPxPyPzE(-p_tot.px(), -p_tot.py(), -p_tot.pz(), p_tot.energy());
    const auto p1_rest = Momentum(pair.at(0)->momentum()).lorentzBoost(p_tot_inv);
    const auto p1_norm = p1_rest.p();
    if (p1_norm <= 0.)
      return false;

    // new momentum in the rest frame (two-body decay kinematics for the new masses)
    const auto lambda = (mass_pair2 - std::pow(m1 + m2, 2)) * (mass_pair2 - std::pow(m1 - m2, 2));
    const auto p_new = std::sqrt(std::max(lambda, 0.)) / (2. * mass_pair);
    const auto scale = p_new / p1_norm;
    const auto p1_new_rest =
        Momentum::fromPxPyPzM(scale * p1_rest.px(), scale * p1_rest.py(), scale * p1_rest.pz(), m1);
    const auto p2_new_rest =
        Momentum::fromPxPyPzM(-scale * p1_rest.px(), -scale * p1_rest.py(), -scale * p1_rest.pz(), m2);
    const auto p1_new = Momentum(p1_new_rest).lorentzBoost(p_tot);  // rest frame -> lab
    const auto p2_new = Momentum(p2_new_rest).lorentzBoost(p_tot);

    // event record consistency: four-momentum of pair must match central system of the record
    auto reference = p_tot;
    if (event.hasRole(Particle::Role::Intermediate))
      reference = event.oneWithRole(Particle::Role::Intermediate).momentum();
    const auto delta_comp = (p1_new + p2_new - reference).pVector(), ref_comp = reference.pVector();
    for (size_t i = 0; i < 4; ++i)
      if (std::fabs(delta_comp[i]) > kMomentumTolerance * (1. + std::fabs(ref_comp[i]))) {
        CG_WARNING("ResonanceMassSmearing:run")
            << "Four-momentum not conserved by the mass smearing (component " << i << ").";
        return false;
      }
    pair.at(0)->setMomentum(p1_new, true);  // particles carry no separate mass
    pair.at(1)->setMomentum(p2_new, true);  // (invariant mass of their four-momentum, kept off-shell here)
    event.metadata["resonanceMassSmearing:m1"] = static_cast<float>(m1);
    event.metadata["resonanceMassSmearing:m2"] = static_cast<float>(m2);
    return true;
  }

private:
  static constexpr double kMomentumTolerance = 1.e-9;  ///< relative four-momentum conservation tolerance
  const std::vector<int> pdg_ids_;
  const double num_widths_;
  const std::unique_ptr<utils::RandomGenerator> random_generator_;
};
REGISTER_MODIFIER("resonanceMassSmearing", ResonanceMassSmearing);
