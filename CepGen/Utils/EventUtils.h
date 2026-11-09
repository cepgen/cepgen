/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#ifndef CepGen_Utils_EventUtils_h
#define CepGen_Utils_EventUtils_h

#include <cmath>

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen::utils {
  /// Generate a standard single-dissociative \f$pp \to p^\ast (\gamma\gamma \to \mu^+\mu^-) p\f$ LPAIR event
  inline Event generateLPAIREvent() {
    auto evt = Event::minimal(2);  // one event with two outgoing particles (leptons)

    // generate positive-z incoming beam kinematics
    auto& ib1 = evt.oneWithRole(Particle::Role::IncomingBeam1);
    ib1.setPdgId(PDG::proton);
    ib1.setMomentum(Momentum::fromPxPyPzE(0., 0., 6.5e3, -1.), false);

    // generate negative-z incoming beam kinematics
    auto& ib2 = evt.oneWithRole(Particle::Role::IncomingBeam2);
    ib2.setPdgId(PDG::proton);
    ib2.setMomentum(Momentum::fromPxPyPzE(0., 0., -6.5e3, -1.), false);

    // generate positive-z outgoing beam kinematics
    auto& ob1 = evt.oneWithRole(Particle::Role::OutgoingBeam1);
    ob1.setPdgId(PDG::proton);
    ob1.setMomentum(Momentum::fromPxPyPzE(-7.875321, 8.186351, 6.403512e3, 6.403704e3), true);

    // generate negative-z outgoing beam kinematics
    auto& ob2 = evt.oneWithRole(Particle::Role::OutgoingBeam2);
    ob2.setPdgId(PDG::proton);
    ob2.setMomentum(Momentum::fromPxPyPzE(-2.725610e-2, 7.565269e-3, -6.425336e3, 6.425336e3), false);

    // generate positive-z incoming photon kinematics
    auto& parton1 = evt.oneWithRole(Particle::Role::Parton1);
    parton1.setPdgId(PDG::photon);
    parton1.setMomentum(Momentum::fromPxPyPzE(7.875321, -8.186351, 9.648800e1, 9.629600e1), true);

    // generate negative-z incoming photon kinematics
    auto& parton2 = evt.oneWithRole(Particle::Role::Parton2);
    parton2.setPdgId(PDG::photon);
    parton2.setMomentum(Momentum::fromPxPyPzE(2.725610e-2, -7.565269e-3, -7.466409e1, 7.466409e1), true);
    evt.oneWithRole(Particle::Role::Intermediate).setMomentum(parton1.momentum() + parton2.momentum(), true);

    // generate two-lepton system kinematics
    auto oc = evt[Particle::Role::CentralSystem];
    oc[0].get().setPdgId(PDG::muon, -1);
    oc[0].get().setMomentum(Momentum::fromPxPyPzE(2.193109e1, -6.725967e1, -4.248568e1, 8.252200e1), false);
    oc[1].get().setPdgId(PDG::muon, +1);
    oc[1].get().setMomentum(Momentum::fromPxPyPzE(-1.402852e1, 5.906575e1, 6.430959e1, 8.843809e1), false);
    return evt;
  }

  /// Generate a four-momentum conserving \f$pp \to p (\gamma\gamma \to X^+X^-) p\f$ event
  /// \param[in] pdg_id Signed PDG identifier of the first central particle (the second one is its anti-particle)
  /// \param[in] cos_theta Cosine of the polar angle of the first central particle in the two-photon rest frame
  /// \param[in] phi Azimuthal angle of the first central particle in the two-photon rest frame
  /// \param[in] e_parton1 Energy of the positive-z incoming photon (in GeV)
  /// \param[in] e_parton2 Energy of the negative-z incoming photon (in GeV)
  inline Event generateBalancedPairEvent(
      spdgid_t pdg_id, double cos_theta = 0.3, double phi = 1., double e_parton1 = 60., double e_parton2 = 40.) {
    constexpr double e_beam = 6.5e3;
    auto evt = Event::minimal(2);
    const auto set = [](Particle& part, spdgid_t id, const Momentum& mom) {
      part.setIntegerPdgId(id);
      part.setMomentum(mom, true);
    };
    const auto pz_beam = std::sqrt(e_beam * e_beam - std::pow(PDG::get().mass(PDG::proton), 2));
    const auto ib1 = Momentum::fromPxPyPzE(0., 0., +pz_beam, e_beam),
               ib2 = Momentum::fromPxPyPzE(0., 0., -pz_beam, e_beam),
               pa1 = Momentum::fromPxPyPzE(0., 0., +e_parton1, e_parton1),
               pa2 = Momentum::fromPxPyPzE(0., 0., -e_parton2, e_parton2);
    set(evt.oneWithRole(Particle::Role::IncomingBeam1), PDG::proton, ib1);
    set(evt.oneWithRole(Particle::Role::IncomingBeam2), PDG::proton, ib2);
    set(evt.oneWithRole(Particle::Role::Parton1), PDG::photon, pa1);
    set(evt.oneWithRole(Particle::Role::Parton2), PDG::photon, pa2);
    set(evt.oneWithRole(Particle::Role::OutgoingBeam1), PDG::proton, ib1 - pa1);
    set(evt.oneWithRole(Particle::Role::OutgoingBeam2), PDG::proton, ib2 - pa2);
    const auto central_system = pa1 + pa2;
    set(evt.oneWithRole(Particle::Role::Intermediate), PDG::photon, central_system);

    // build the pair back-to-back in the two-photon rest frame, then boost it along the z-axis
    const auto mass = central_system.mass(), gamma = central_system.energy() / mass,
               beta_gamma = central_system.pz() / mass, energy_cm = 0.5 * mass,
               p_cm = std::sqrt(energy_cm * energy_cm - std::pow(PDG::get().mass(pdg_id), 2)),
               sin_theta = std::sqrt(1. - cos_theta * cos_theta);
    auto oc = evt[Particle::Role::CentralSystem];
    for (size_t i = 0; i < 2; ++i) {
      const double sign = i == 0 ? +1. : -1., pz_cm = sign * p_cm * cos_theta;
      set(oc[i].get(),
          i == 0 ? pdg_id : -pdg_id,
          Momentum::fromPxPyPzE(sign * p_cm * sin_theta * std::cos(phi),
                                sign * p_cm * sin_theta * std::sin(phi),
                                gamma * pz_cm + beta_gamma * energy_cm,
                                gamma * energy_cm + beta_gamma * pz_cm));
    }
    return evt;
  }
}  // namespace cepgen::utils

#endif
