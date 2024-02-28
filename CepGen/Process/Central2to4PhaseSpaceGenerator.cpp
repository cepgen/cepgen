/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2024  Laurent Forthomme
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

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Central2to4PhaseSpaceGenerator.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  Central2to4PhaseSpaceGenerator::Central2to4PhaseSpaceGenerator(const ParametersList& params)
      : CentralPhaseSpaceGenerator(params) {}

  void Central2to4PhaseSpaceGenerator::initialise() {
    const auto& kin_cuts = process().kinematics().cuts().central;
    const auto lim_rap = kin_cuts.rapidity_single.truncate(Limits{-6., 6.});
    process()
        .defineVariable(m_y_c1_, proc::Process::Mapping::linear, lim_rap, "y1", "First outgoing particle rapidity")
        .defineVariable(m_y_c2_, proc::Process::Mapping::linear, lim_rap, "y2", "Second outgoing particle rapidity")
        .defineVariable(m_pt_diff_,
                        proc::Process::Mapping::linear,
                        kin_cuts.pt_diff.truncate(Limits{0., 500.}),
                        "pt_diff",
                        "Final state particles transverse momentum difference")
        .defineVariable(m_phi_pt_diff_,
                        proc::Process::Mapping::linear,
                        kin_cuts.phi_diff.truncate(Limits{0., 2. * M_PI}),
                        "phi_pt_diff",
                        "Final state particles azimuthal angle difference");
  }

  double Central2to4PhaseSpaceGenerator::generateKinematics() {
    const auto& kin = process().kinematics();
    if (!kin.cuts().central.rapidity_diff.contains(std::fabs(m_y_c1_ - m_y_c2_)))  // rapidity distance
      return 0.;
    {
      const auto qt_sum = (process().q1() + process().q2()).transverse();  // two-parton system
      const auto pt_diff = Momentum::fromPtEtaPhiE(m_pt_diff_, 0., m_phi_pt_diff_);
      const auto pt_c1 = 0.5 * (qt_sum + pt_diff), pt_c2 = 0.5 * (qt_sum - pt_diff);
      const auto p1t = pt_c1.pt(), p2t = pt_c2.pt();
      // apply user cuts on central system
      if (!kin.cuts().central.pt_single.contains(p1t) || !single_limits_.pt_single.contains(p1t))
        return 0.;
      if (!kin.cuts().central.pt_single.contains(p2t) || !single_limits_.pt_single.contains(p2t))
        return 0.;
      if (!kin.cuts().central.pt_diff.contains(std::fabs(p1t - p2t)))  // transverse momentum difference
        return 0.;
      //--- four-momenta of the outgoing central particles
      process().pc(0) = Momentum::fromPtYPhiM(p1t, m_y_c1_, pt_c1.phi(), PDG::get().mass(particles_.at(0)));
      process().pc(1) = Momentum::fromPtYPhiM(p2t, m_y_c2_, pt_c2.phi(), PDG::get().mass(particles_.at(1)));
    }

    //--- window in central system invariant mass
    const auto invm = (process().pc(0) + process().pc(1)).mass();
    if (!kin.cuts().central.mass_sum.contains(invm))
      return 0.;

    //--- compute and sanitise the momentum losses
    const auto amt1 = process().pc(0).massT() / process().sqrtS(), amt2 = process().pc(1).massT() / process().sqrtS();
    static const auto x_lim = Limits{0., 1.};
    const auto x1 = amt1 * exp(+m_y_c1_) + amt2 * exp(+m_y_c2_);
    if (!x_lim.contains(x1))
      return 0.;
    const auto x2 = amt1 * exp(-m_y_c1_) + amt2 * exp(-m_y_c2_);
    if (!x_lim.contains(x2))
      return 0.;

    //--- additional conditions for energy-momentum conservation
    const auto s = process().s(), mx2 = process().mX2(), my2 = process().mY2();
    if (!kin.incomingBeams().positive().elastic() && x2 * s - invm - process().q2().p2() <= mx2)
      return 0.;
    if (!kin.incomingBeams().negative().elastic() && x1 * s - invm - process().q1().p2() <= my2)
      return 0.;

    //--- four-momenta of the outgoing protons (or remnants)

    const auto px_p = (1. - x1) * process().pA().p() * M_SQRT2, px_m = (mx2 + process().q1().p2()) * 0.5 / px_p;
    const auto py_m = (1. - x2) * process().pB().p() * M_SQRT2, py_p = (my2 + process().q2().p2()) * 0.5 / py_m;
    CG_DEBUG_LOOP("2to4:pxy") << "px+ = " << px_p << " / px- = " << px_m << "\n\t"
                              << "py+ = " << py_p << " / py- = " << py_m << ".";

    const auto px = -Momentum(process().q1()).setPz((px_p - px_m) * M_SQRT1_2).setEnergy((px_p + px_m) * M_SQRT1_2),
               py = -Momentum(process().q2()).setPz((py_p - py_m) * M_SQRT1_2).setEnergy((py_p + py_m) * M_SQRT1_2);

    CG_DEBUG_LOOP("2to4:remnants") << "First remnant:  " << px << ", mass = " << px.mass() << "\n\t"
                                   << "Second remnant: " << py << ", mass = " << py.mass() << ".";

    if (std::fabs(px.mass2() - mx2) > NUM_LIMITS) {
      CG_WARNING("2to4:px") << "Invalid X system squared mass: " << px.mass2() << "/" << mx2 << ".";
      return 0.;
    }
    if (std::fabs(py.mass2() - my2) > NUM_LIMITS) {
      CG_WARNING("2to4:py") << "Invalid Y system squared mass: " << py.mass2() << "/" << my2 << ".";
      return 0.;
    }

    //--- four-momenta of the intermediate partons
    const double norm = 1. / process().wCM() / process().wCM() / s, prefac = 0.5 / std::sqrt(norm);
    {  // positive-z incoming parton collinear kinematics
      const double tau1 = norm * process().q1().p2() / x1 / x1;
      process().q1().setPz(+prefac * x1 * (1. - tau1)).setEnergy(+prefac * x1 * (1. + tau1));
    }
    {  // negative-z incoming parton collinear kinematics
      const double tau2 = norm * process().q2().p2() / x2 / x2;
      process().q2().setPz(-prefac * x2 * (1. - tau2)).setEnergy(+prefac * x2 * (1. + tau2));
    }

    CG_DEBUG_LOOP("2to4:partons") << "Squared c.m. energy = " << s << " GeV^2\n\t"
                                  << "First parton: " << process().q1() << ", mass2 = " << process().q1().mass2()
                                  << ", x1 = " << x1 << ", p = " << process().q1().p() << "\n\t"
                                  << "Second parton: " << process().q2() << ", mass2 = " << process().q2().mass2()
                                  << ", x2 = " << x2 << ", p = " << process().q2().p() << ".";

    // randomise the charge of outgoing system
    const short sign = process().randomGenerator().uniformInt(0, 1) == 1 ? 1 : -1;
    process().event()[Particle::CentralSystem][0].get().setChargeSign(+sign).setStatus(Particle::Status::FinalState);
    process().event()[Particle::CentralSystem][1].get().setChargeSign(-sign).setStatus(Particle::Status::FinalState);
    process().x1() = x1;
    process().x2() = x2;
    process().pX() = px;
    process().pY() = py;
    return prefactor_ * m_pt_diff_;
  }

  ParametersDescription Central2to4PhaseSpaceGenerator::description() {
    auto desc = CentralPhaseSpaceGenerator::description();
    desc.setDescription("2-to-4 process");
    return desc;
  }
}  // namespace cepgen
