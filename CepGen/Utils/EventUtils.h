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
}  // namespace cepgen::utils

#endif
