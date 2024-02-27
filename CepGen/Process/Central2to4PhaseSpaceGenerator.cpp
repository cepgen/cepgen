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
#include "CepGen/Modules/CentralPhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/CentralPhaseSpaceGenerator.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  /// A 2-to-4 (or 2-to-2 central) phase space generator
  class Central2to4PhaseSpaceGenerator : public CentralPhaseSpaceGenerator {
  public:
    explicit Central2to4PhaseSpaceGenerator(const ParametersList& params)
        : CentralPhaseSpaceGenerator(params), particles_(steer<pdgids_t>("ids")) {}

    size_t ndim() const override { return 4; }
    const pdgids_t& particles() const override { return particles_; }

    void initialise() override {
      const auto& kin_cuts = process().kinematics().cuts().central;
      const auto lim_rap = kin_cuts.rapidity_single.truncate(Limits{-6., 6.});
      process().defineVariable(
          m_y_c1_, proc::Process::Mapping::linear, lim_rap, "y1", "First outgoing particle rapidity");
      process().defineVariable(
          m_y_c2_, proc::Process::Mapping::linear, lim_rap, "y2", "Second outgoing particle rapidity");

      const auto lim_pt_diff = kin_cuts.pt_diff.truncate(Limits{0., 500.});
      process().defineVariable(m_pt_diff_,
                               proc::Process::Mapping::linear,
                               lim_pt_diff,
                               "pt_diff",
                               "Final state particles transverse momentum difference");

      const auto lim_phi_diff = kin_cuts.phi_diff.truncate(Limits{0., 2. * M_PI});
      process().defineVariable(m_phi_pt_diff_,
                               proc::Process::Mapping::linear,
                               lim_phi_diff,
                               "phi_pt_diff",
                               "Final state particles azimuthal angle difference");
    }

    double generateKinematics() override {
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
      if (process().x1() = amt1 * exp(+m_y_c1_) + amt2 * exp(+m_y_c2_); !x_lim.contains(process().x1()))
        return 0.;
      if (process().x2() = amt1 * exp(-m_y_c1_) + amt2 * exp(-m_y_c2_); !x_lim.contains(process().x2()))
        return 0.;

      //--- additional conditions for energy-momentum conservation
      if (!kin.incomingBeams().positive().elastic() &&
          std::sqrt(process().x2() * process().s() - invm - process().q2().p2()) <= process().mX())
        return 0.;
      if (!kin.incomingBeams().negative().elastic() &&
          std::sqrt(process().x1() * process().s() - invm - process().q1().p2()) <= process().mY())
        return 0.;

      //--- four-momenta of the outgoing protons (or remnants)

      const auto px_p = (1. - process().x1()) * process().pA().p() * M_SQRT2,
                 px_m = (process().mX2() + process().q1().p2()) * 0.5 / px_p;
      const auto py_m = (1. - process().x2()) * process().pB().p() * M_SQRT2,
                 py_p = (process().mY2() + process().q2().p2()) * 0.5 / py_m;
      CG_DEBUG_LOOP("2to4:pxy") << "px+ = " << px_p << " / px- = " << px_m << "\n\t"
                                << "py+ = " << py_p << " / py- = " << py_m << ".";

      process().pX() = -Momentum(process().q1()).setPz((px_p - px_m) * M_SQRT1_2).setEnergy((px_p + px_m) * M_SQRT1_2);
      process().pY() = -Momentum(process().q2()).setPz((py_p - py_m) * M_SQRT1_2).setEnergy((py_p + py_m) * M_SQRT1_2);

      CG_DEBUG_LOOP("2to4:remnants") << "First remnant:  " << process().pX() << ", mass = " << process().pX().mass()
                                     << "\n\t"
                                     << "Second remnant: " << process().pY() << ", mass = " << process().pY().mass()
                                     << ".";

      if (std::fabs(process().pX().mass2() - process().mX2()) > NUM_LIMITS) {
        CG_WARNING("2to4:px") << "Invalid X system squared mass: " << process().pX().mass2() << "/" << process().mX2()
                              << ".";
        return 0.;
      }
      if (std::fabs(process().pY().mass2() - process().mY2()) > NUM_LIMITS) {
        CG_WARNING("2to4:py") << "Invalid Y system squared mass: " << process().pY().mass2() << "/" << process().mY2()
                              << ".";
        return 0.;
      }

      //--- four-momenta of the intermediate partons
      const double norm = 1. / process().wCM() / process().wCM() / process().s(),
                   prefac = 0.5 * process().wCM() * process().sqrtS();
      {  // positive-z incoming parton collinear kinematics
        const double tau1 = norm * process().q1().p2() / process().x1() / process().x1();
        process().q1().setPz(+prefac * process().x1() * (1. - tau1)).setEnergy(+prefac * process().x1() * (1. + tau1));
      }
      {  // negative-z incoming parton collinear kinematics
        const double tau2 = norm * process().q2().p2() / process().x2() / process().x2();
        process().q2().setPz(-prefac * process().x2() * (1. - tau2)).setEnergy(+prefac * process().x2() * (1. + tau2));
      }

      CG_DEBUG_LOOP("2to4:partons") << "Squared c.m. energy = " << process().s() << " GeV^2\n\t"
                                    << "First parton: " << process().q1() << ", mass2 = " << process().q1().mass2()
                                    << ", x1 = " << process().x1() << ", p = " << process().q1().p() << "\n\t"
                                    << "Second parton: " << process().q2() << ", mass2 = " << process().q2().mass2()
                                    << ", x2 = " << process().x2() << ", p = " << process().q2().p() << ".";

      // randomise the charge of outgoing system
      const short sign = process().randomGenerator().uniformInt(0, 1) == 1 ? 1 : -1;
      process().event()[Particle::CentralSystem][0].get().setChargeSign(+sign).setStatus(Particle::Status::FinalState);
      process().event()[Particle::CentralSystem][1].get().setChargeSign(-sign).setStatus(Particle::Status::FinalState);
      return prefactor_ * m_pt_diff_;
    }

  protected:
    // mapped variables
    double m_y_c1_{0.};         ///< Rapidity of the first central particle
    double m_y_c2_{0.};         ///< Rapidity of the second central particle
    double m_pt_diff_{0.};      ///< Transverse momentum difference for the two central particle
    double m_phi_pt_diff_{0.};  ///< Azimuthal angle difference for the two central particles

  private:
    // factor 1/4 from jacobian of transformations
    static constexpr double prefactor_ = 0.25 * 0.0625 * M_1_PI * M_1_PI;
    const pdgids_t particles_;  ///< Type of particles produced in the final state
  };
}  // namespace cepgen
REGISTER_CENTRAL_PSGEN("2to4", Central2to4PhaseSpaceGenerator);
