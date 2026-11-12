/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2026  Laurent Forthomme
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
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace cepgen::hadr;

Hadroniser::Hadroniser(const ParametersList& plist)
    : EventModifier(plist), fragment_remnants_(steer<bool>("remnantsFragmentation")) {}

ParametersDescription Hadroniser::description() {
  auto desc = EventModifier::description();
  desc.add("remnantsFragmentation", true).setDescription("Apply the fragmentation algorithm to proton remnants");
  return desc;
}

Hadroniser::DiffractiveState::DiffractiveState(utils::RandomGenerator& rng,
                                               long beam_pdg_id,
                                               const PartonsContent& partons_content)
    : random_generator_(&rng), beam_pdg_id_(beam_pdg_id) {
  if (std::labs(beam_pdg_id_) != PDG::proton)
    throw CG_FATAL("Hadroniser:DiffractiveState")
        << "Beam particle with PDG id " << beam_pdg_id_
        << " is not supported. Only protons and antiprotons can be dissociated.";
  if (partons_content.empty())
    throw CG_FATAL("Hadroniser:DiffractiveState") << "Flavour content must be provided.";

  // shoot valence flavour content
  double total_weight = 0.;
  for (const auto& content : partons_content)
    total_weight += content.weight;
  const long sign = beam_pdg_id_ < 0 ? -1 : +1;
  auto draw = random_generator_->uniform(0., total_weight);
  for (const auto& content : partons_content) {
    valquark_pdg_id_ = sign * content.valquark_pdg_id;
    diquark_pdg_id_ = sign * content.diquark_pdg_id;
    if ((draw -= content.weight) < 0.)
      break;
  }
}

void Hadroniser::DiffractiveState::initialise(const Momentum& system_momentum,
                                              const Momentum& beam_momentum,
                                              double primordial_kt) {
  system_mass_ = system_momentum.mass();
  if (valquark_mass_ <= 0.)
    valquark_mass_ = massFromPdgId(valquark_pdg_id_);
  if (diquark_mass_ <= 0.)
    diquark_mass_ = massFromPdgId(diquark_pdg_id_);

  const auto mass2 = system_mass_ * system_mass_;
  const auto sum = valquark_mass_ + diquark_mass_, diff = valquark_mass_ - diquark_mass_;
  const auto lambda = (mass2 - sum * sum) * (mass2 - diff * diff);
  if (lambda <= 0. || system_mass_ <= 0.) {
    CG_DEBUG("pythia8:DiffractiveState") << "Diffractive system mass unphysical. System invalid for fragmentation.";
    valid_ = false;
    return;
  }
  const auto p_star = 0.5 * std::sqrt(lambda) / system_mass_;

  // the spectator follows the beam direction, as seen from the rest frame of the dissociated system
  const auto beam_in_system =
      Momentum(beam_momentum)
          .lorentzBoost(Momentum::fromPxPyPzE(
              -system_momentum.px(), -system_momentum.py(), -system_momentum.pz(), system_momentum.energy()));
  const auto axis = beam_in_system.p() > 0.
                        ? beam_in_system.unit()
                        : Momentum::randomThetaPhi(*random_generator_, {-M_PI, M_PI}, {0., 2. * M_PI});

  auto direction = axis;
  if (primordial_kt > 0.) {
    const auto px = random_generator_->normal(0., primordial_kt * M_SQRT1_2),
               py = random_generator_->normal(0., primordial_kt * M_SQRT1_2);
    const auto pt = std::hypot(px, py);
    if (pt > 0.) {
      const auto sin_theta = std::min(pt / p_star, 1.), cos_theta = std::sqrt(1. - sin_theta * sin_theta);
      const auto helper =
          std::fabs(axis.px()) < 0.9 ? Momentum::fromPxPyPzE(1., 0., 0., 1.) : Momentum::fromPxPyPzE(0., 1., 0., 1.);
      const auto e1 = axis.crossProduct(helper).unit(), e2 = axis.crossProduct(e1);
      direction = cos_theta * axis + sin_theta * (px / pt * e1 + py / pt * e2);
    }
  }
  diquark_momentum_ = p_star * direction;
}

bool Hadroniser::DiffractiveState::momentumConserved(const Fragments& fragments, double tolerance) const {
  Momentum sum{0., 0., 0., system_mass_};
  for (const auto& fragment : fragments)
    if (fragment.stable)
      sum -= fragment.momentum;
  if (const auto residual = sum.p(); residual > tolerance * system_mass_) {
    CG_WARNING("DiffractiveState:momentumConserved")
        << "Momentum not conserved down to " << tolerance << " tolerance after hadronisation (residual=" << residual
        << " GeV for a system of mass " << system_mass_ << " GeV).";
    return false;
  }
  return true;
}

double Hadroniser::DiffractiveState::massFromPdgId(long pdg_id) const { return PDG::get().mass(pdg_id); }
