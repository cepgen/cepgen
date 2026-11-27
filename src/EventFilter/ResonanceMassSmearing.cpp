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

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Message.h"
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
    for (auto& part_ref : event[Particle::Role::CentralSystem]) {
      auto& part = part_ref.get();
      if (std::find(pdg_ids_.begin(), pdg_ids_.end(), std::abs(static_cast<int>(part.integerPdgId()))) !=
          pdg_ids_.end())
        pair.emplace_back(&part);
    }
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

    // pair four-momentum and invariant mass
    std::array<double, 4> p_tot{};  // (px, py, pz, E)
    for (const auto* part : pair) {
      const auto& mom = part->momentum();
      p_tot[0] += mom.px(), p_tot[1] += mom.py(), p_tot[2] += mom.pz(), p_tot[3] += mom.energy();
    }
    const auto mass_pair2 = p_tot[3] * p_tot[3] - p_tot[0] * p_tot[0] - p_tot[1] * p_tot[1] - p_tot[2] * p_tot[2];
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

    // pair rest frame: boost velocity
    const std::array<double, 3> beta{p_tot[0] / p_tot[3], p_tot[1] / p_tot[3], p_tot[2] / p_tot[3]};
    const auto beta2 = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
    const auto gamma = 1. / std::sqrt(1. - beta2);
    const auto boost = [&](const std::array<double, 4>& p, double sign) {  // sign = +1: lab -> rest, -1: rest -> lab
      const auto bp = beta[0] * p[0] + beta[1] * p[1] + beta[2] * p[2];
      const auto gamma2 = beta2 > 0. ? (gamma - 1.) / beta2 : 0.;
      const auto coef = sign > 0 ? gamma2 * bp - gamma * p[3] : gamma2 * bp + gamma * p[3];
      return std::array<double, 4>{
          p[0] + coef * beta[0], p[1] + coef * beta[1], p[2] + coef * beta[2], gamma * (p[3] - sign * bp)};
    };

    const auto& mom1 = pair.at(0)->momentum();
    const auto p1_rest = boost({mom1.px(), mom1.py(), mom1.pz(), mom1.energy()}, +1.);
    const auto p1_norm = std::sqrt(p1_rest[0] * p1_rest[0] + p1_rest[1] * p1_rest[1] + p1_rest[2] * p1_rest[2]);
    if (p1_norm <= 0.)
      return false;

    // new momentum in the rest frame (two-body decay kinematics for the new masses)
    const auto lambda = (mass_pair2 - std::pow(m1 + m2, 2)) * (mass_pair2 - std::pow(m1 - m2, 2));
    const auto p_new = std::sqrt(std::max(lambda, 0.)) / (2. * mass_pair);
    const std::array<double, 4> p1_new_rest{p_new * p1_rest[0] / p1_norm,
                                            p_new * p1_rest[1] / p1_norm,
                                            p_new * p1_rest[2] / p1_norm,
                                            std::sqrt(p_new * p_new + m1 * m1)},
        p2_new_rest{-p1_new_rest[0], -p1_new_rest[1], -p1_new_rest[2], std::sqrt(p_new * p_new + m2 * m2)};
    const auto p1_new = boost(p1_new_rest, -1.), p2_new = boost(p2_new_rest, -1.);

    // event record consistency: the pair four-momentum (hence the intermediate two-parton system, the beam partons and
    // the outgoing beams) must be preserved, to the precision the event record parentage check relies on
    for (size_t i = 0; i < 4; ++i)
      if (std::fabs(p1_new[i] + p2_new[i] - p_tot[i]) > kMomentumTolerance * (1. + std::fabs(p_tot[i]))) {
        CG_WARNING("ResonanceMassSmearing:run")
            << "Four-momentum not conserved by the mass smearing (component " << i << ").";
        return false;
      }
    // the particles carry no separate mass: it is the invariant mass of their four-momentum, kept off-shell here
    pair.at(0)->setMomentum(Momentum::fromPxPyPzE(p1_new[0], p1_new[1], p1_new[2], p1_new[3]), true);
    pair.at(1)->setMomentum(Momentum::fromPxPyPzE(p2_new[0], p2_new[1], p2_new[2], p2_new[3]), true);
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
