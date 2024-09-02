/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  PhaseSpaceGenerator::PhaseSpaceGenerator(const ParametersList& params) : NamedModule(params) {}

  bool PhaseSpaceGenerator::constrainBeamKinematics(proc::FactorisedProcess* proc) {
    const auto& kin = proc->kinematics();
    //--- window in central system invariant mass
    const auto invariant_mass = (proc->q1() + proc->q2()).mass();
    if (!kin.cuts().central.mass_sum.contains(invariant_mass))
      return 0.;

    //--- compute and sanitise the momentum losses
    static const auto x_lim = Limits{0., 1.};
    proc->x1() = proc->x2() = 0.;
    for (size_t i = 0; i < central().size(); ++i) {
      const auto amt = proc->pc(i).massT() * proc->inverseSqrtS(), ay = proc->pc(i).rapidity();
      proc->x1() += amt * std::exp(+ay);
      proc->x2() += amt * std::exp(-ay);
    }
    if (!x_lim.contains(proc->x1()) || !x_lim.contains(proc->x2()))
      return 0.;
    //--- additional conditions for energy-momentum conservation
    if (!kin.incomingBeams().positive().elastic() &&
        proc->x2() * proc->s() - invariant_mass - proc->q2().p2() <= proc->mX2())
      return 0.;
    if (!kin.incomingBeams().negative().elastic() &&
        proc->x1() * proc->s() - invariant_mass - proc->q1().p2() <= proc->mY2())
      return 0.;

    //--- four-momenta of the outgoing protons (or remnants)

    const auto px_p = (1. - proc->x1()) * proc->pA().p() * M_SQRT2, px_m = (proc->mX2() + proc->q1().p2()) * 0.5 / px_p;
    const auto py_m = (1. - proc->x2()) * proc->pB().p() * M_SQRT2, py_p = (proc->mY2() + proc->q2().p2()) * 0.5 / py_m;
    CG_DEBUG_LOOP("PhaseSpaceGenerator:pxy") << "px+ = " << px_p << " / px- = " << px_m << "\n\t"
                                             << "py+ = " << py_p << " / py- = " << py_m << ".";

    proc->pX() = -Momentum(proc->q1()).setPz((px_p - px_m) * M_SQRT1_2).setEnergy((px_p + px_m) * M_SQRT1_2);
    proc->pY() = -Momentum(proc->q2()).setPz((py_p - py_m) * M_SQRT1_2).setEnergy((py_p + py_m) * M_SQRT1_2);

    CG_DEBUG_LOOP("PhaseSpaceGenerator:remnants")
        << "First remnant:  " << proc->pX() << ", mass = " << proc->pX().mass() << "\n\t"
        << "Second remnant: " << proc->pY() << ", mass = " << proc->pY().mass() << ".";

    if (std::fabs(proc->pX().mass2() - proc->mX2()) > NUM_LIMITS) {
      CG_WARNING("PhaseSpaceGenerator:px")
          << "Invalid X system squared mass: " << proc->pX().mass2() << "/" << proc->mX2() << ".";
      return false;
    }
    if (std::fabs(proc->pY().mass2() - proc->mY2()) > NUM_LIMITS) {
      CG_WARNING("PhaseSpaceGenerator:py")
          << "Invalid Y system squared mass: " << proc->pY().mass2() << "/" << proc->mY2() << ".";
      return false;
    }

    //--- four-momenta of the intermediate partons
    const double norm = 1. / proc->wCM() / proc->wCM() * proc->inverseS(), prefactor = 0.5 / std::sqrt(norm);
    {  // positive-z incoming parton collinear kinematics
      const double tau1 = norm * proc->q1().p2() / proc->x1();
      proc->q1().setPz(+prefactor * (proc->x1() - tau1)).setEnergy(+prefactor * (proc->x1() + tau1));
    }
    {  // negative-z incoming parton collinear kinematics
      const double tau2 = norm * proc->q2().p2() / proc->x2();
      proc->q2().setPz(-prefactor * (proc->x2() - tau2)).setEnergy(+prefactor * (proc->x2() + tau2));
    }

    CG_DEBUG_LOOP("PhaseSpaceGenerator:partons")
        << "Squared c.m. energy = " << proc->s() << " GeV^2\n\t"
        << "First parton: " << proc->q1() << ", mass2 = " << proc->q1().mass2() << ", x1 = " << proc->x1()
        << ", p = " << proc->q1().p() << "\n\t"
        << "Second parton: " << proc->q2() << ", mass2 = " << proc->q2().mass2() << ", x2 = " << proc->x2()
        << ", p = " << proc->q2().p() << ".";

    return true;
  }
}  // namespace cepgen
