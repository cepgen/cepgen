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

#include <memory>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Process/PartonsCollinearPhaseSpaceGenerator.h"
#include "CepGen/Process/PartonsKTPhaseSpaceGenerator.h"
#include "CepGen/Process/PartonsPhaseSpaceGenerator.h"
#include "CepGen/Process/PhaseSpaceGenerator.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  /// A 2-to-4 (or 2-to-2 central) phase space generator
  template <typename T>
  class PhaseSpaceGenerator2to4 : public PhaseSpaceGenerator {
  public:
    explicit PhaseSpaceGenerator2to4(const ParametersList& params)
        : PhaseSpaceGenerator(params),
          part_psgen_(new T(params)),
          int_particles_(steer<std::vector<int> >("ids")),
          particles_(int_particles_.begin(), int_particles_.end()) {}

    static ParametersDescription description() {
      auto desc = PhaseSpaceGenerator::description();
      desc.setDescription("2-to-4 phase space mapper (" + T::description().description() + "/" +
                          T::description().description() + ")");
      desc.add<std::vector<int> >("ids", {}).setDescription("list of particles produced");
      desc += T::description();
      return desc;
    }

    bool ktFactorised() const override {
      CG_ASSERT(part_psgen_);
      return part_psgen_->ktFactorised();
    }

    void setCentralCuts(const cuts::Central& single) override { single_limits_ = single; }

    void initialise(proc::FactorisedProcess* process) override {
      proc_ = process;
      CG_ASSERT(part_psgen_);
      part_psgen_->initialise(process);
      const auto& kin_cuts = proc_->kinematics().cuts().central;
      const auto lim_rap = kin_cuts.rapidity_single.truncate(Limits{-6., 6.});
      (*proc_)
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

    double generate() override {
      CG_ASSERT(part_psgen_);
      if (!part_psgen_->generatePartonKinematics())
        return 0.;
      const auto cent_weight = generateCentralKinematics();
      if (!utils::positive(cent_weight))
        return 0.;
      const auto fluxes_weight = part_psgen_->fluxes();
      if (!utils::positive(fluxes_weight))
        return 0.;
      return fluxes_weight * cent_weight;
    }

    pdgids_t partons() const override {
      CG_ASSERT(part_psgen_);
      return pdgids_t{part_psgen_->positiveFlux().partonPdgId(), part_psgen_->negativeFlux().partonPdgId()};
    }
    pdgids_t central() const override { return particles_; }

  private:
    double generateCentralKinematics() {
      const auto& kin = proc_->kinematics();
      if (!kin.cuts().central.rapidity_diff.contains(std::fabs(m_y_c1_ - m_y_c2_)))  // rapidity distance
        return 0.;
      {
        const auto qt_sum = (proc_->q1() + proc_->q2()).transverse();  // two-parton system
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
        proc_->pc(0) = Momentum::fromPtYPhiM(p1t, m_y_c1_, pt_c1.phi(), PDG::get().mass(particles_.at(0)));
        proc_->pc(1) = Momentum::fromPtYPhiM(p2t, m_y_c2_, pt_c2.phi(), PDG::get().mass(particles_.at(1)));
      }

      //--- window in central system invariant mass
      const auto invm = (proc_->pc(0) + proc_->pc(1)).mass();
      if (!kin.cuts().central.mass_sum.contains(invm))
        return 0.;

      //--- compute and sanitise the momentum losses
      const auto amt1 = proc_->pc(0).massT() / proc_->sqrtS(), amt2 = proc_->pc(1).massT() / proc_->sqrtS();
      static const auto x_lim = Limits{0., 1.};
      const auto x1 = amt1 * exp(+m_y_c1_) + amt2 * exp(+m_y_c2_);
      if (!x_lim.contains(x1))
        return 0.;
      const auto x2 = amt1 * exp(-m_y_c1_) + amt2 * exp(-m_y_c2_);
      if (!x_lim.contains(x2))
        return 0.;

      //--- additional conditions for energy-momentum conservation
      const auto s = proc_->s(), mx2 = proc_->mX2(), my2 = proc_->mY2();
      if (!kin.incomingBeams().positive().elastic() && x2 * s - invm - proc_->q2().p2() <= mx2)
        return 0.;
      if (!kin.incomingBeams().negative().elastic() && x1 * s - invm - proc_->q1().p2() <= my2)
        return 0.;

      //--- four-momenta of the outgoing protons (or remnants)

      const auto px_p = (1. - x1) * proc_->pA().p() * M_SQRT2, px_m = (mx2 + proc_->q1().p2()) * 0.5 / px_p;
      const auto py_m = (1. - x2) * proc_->pB().p() * M_SQRT2, py_p = (my2 + proc_->q2().p2()) * 0.5 / py_m;
      CG_DEBUG_LOOP("2to4:pxy") << "px+ = " << px_p << " / px- = " << px_m << "\n\t"
                                << "py+ = " << py_p << " / py- = " << py_m << ".";

      const auto px = -Momentum(proc_->q1()).setPz((px_p - px_m) * M_SQRT1_2).setEnergy((px_p + px_m) * M_SQRT1_2),
                 py = -Momentum(proc_->q2()).setPz((py_p - py_m) * M_SQRT1_2).setEnergy((py_p + py_m) * M_SQRT1_2);

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
      const double norm = 1. / proc_->wCM() / proc_->wCM() / s, prefac = 0.5 / std::sqrt(norm);
      {  // positive-z incoming parton collinear kinematics
        const double tau1 = norm * proc_->q1().p2() / x1 / x1;
        proc_->q1().setPz(+prefac * x1 * (1. - tau1)).setEnergy(+prefac * x1 * (1. + tau1));
      }
      {  // negative-z incoming parton collinear kinematics
        const double tau2 = norm * proc_->q2().p2() / x2 / x2;
        proc_->q2().setPz(-prefac * x2 * (1. - tau2)).setEnergy(+prefac * x2 * (1. + tau2));
      }

      CG_DEBUG_LOOP("2to4:partons") << "Squared c.m. energy = " << s << " GeV^2\n\t"
                                    << "First parton: " << proc_->q1() << ", mass2 = " << proc_->q1().mass2()
                                    << ", x1 = " << x1 << ", p = " << proc_->q1().p() << "\n\t"
                                    << "Second parton: " << proc_->q2() << ", mass2 = " << proc_->q2().mass2()
                                    << ", x2 = " << x2 << ", p = " << proc_->q2().p() << ".";

      // randomise the charge of outgoing system
      const short sign = proc_->randomGenerator().uniformInt(0, 1) == 1 ? 1 : -1;
      proc_->event()[Particle::CentralSystem][0].get().setChargeSign(+sign).setStatus(Particle::Status::FinalState);
      proc_->event()[Particle::CentralSystem][1].get().setChargeSign(-sign).setStatus(Particle::Status::FinalState);
      proc_->x1() = x1;
      proc_->x2() = x2;
      proc_->pX() = px;
      proc_->pY() = py;
      return prefactor_ * m_pt_diff_;
    }

    // factor 1/4 from jacobian of transformations
    static constexpr double prefactor_ = 0.25 * 0.0625 * M_1_PI * M_1_PI;
    static constexpr double NUM_LIMITS = 1.e-3;  ///< Numerical limits for sanity comparisons (MeV/mm-level)

    const std::unique_ptr<PartonsPhaseSpaceGenerator> part_psgen_;
    const std::vector<int> int_particles_;  ///< Type of particles produced in the final state (integer values)
    const pdgids_t particles_;              ///< Type of particles produced in the final state (PDG ids)

    proc::FactorisedProcess* proc_{nullptr};  //NOT owning

    cuts::Central single_limits_;  ///< Limits to be applied on single central system's particles
    // mapped variables
    double m_y_c1_{0.};         ///< Rapidity of the first central particle
    double m_y_c2_{0.};         ///< Rapidity of the second central particle
    double m_pt_diff_{0.};      ///< Transverse momentum difference for the two central particle
    double m_phi_pt_diff_{0.};  ///< Azimuthal angle difference for the two central particles
  };

  typedef PhaseSpaceGenerator2to4<PartonsKTPhaseSpaceGenerator> KT2to4;
  typedef PhaseSpaceGenerator2to4<PartonsCollinearPhaseSpaceGenerator> Coll2to4;
}  // namespace cepgen
REGISTER_PSGEN("kt2to4", KT2to4);
REGISTER_PSGEN("coll2to4", Coll2to4);
