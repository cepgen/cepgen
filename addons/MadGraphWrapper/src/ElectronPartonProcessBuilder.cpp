/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025-2026  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGenMadGraph/Process.h"
#include "CepGenMadGraph/ProcessBuilder.h"

namespace cepgen::mg5amc {
  class ElectronPartonProcessBuilder final : public ProcessBuilder {
  public:
    explicit ElectronPartonProcessBuilder(const ParametersList& params, bool load_library = true)
        : ProcessBuilder(params, load_library) {
      auto central_system = process().centralSystem();
      if (std::abs(*central_system.begin()) == PDG::electron) {
        electron_first_ = true;
        central_system.erase(central_system.begin());  // first particle is lepton
      } else if (std::abs(*central_system.rbegin()) == PDG::electron) {
        electron_first_ = false;
        central_system.pop_back();  // last particle is lepton
      } else
        throw CG_FATAL("ElectronPartonProcessBuilder")
            << "No electron/positron found in mg5_aMC process particles list.";
      setCentral(central_system);  // electron/positron stripped off central system
    }

    proc::ProcessPtr clone() const override {
      return std::make_unique<ElectronPartonProcessBuilder>(parameters(), false);
    }

    static ParametersDescription description() {
      auto desc = ProcessBuilder::description();
      desc.setDescription("MadGraph_aMC electron-parton process builder");
      return desc;
    }

    bool fluxlessBeam(PartonSampler::Side side) const override {
      return side == (electron_first_ ? PartonSampler::Side::positive : PartonSampler::Side::negative);
    }

    void prepareFactorisedPhaseSpace() override {
      // the hadron-side parton is the last (first) incoming MadGraph parton if the lepton is first (last)
      const auto& mg_partons = process().intermediatePartons();
      const auto mg_hadron_parton = electron_first_ ? *mg_partons.rbegin() : *mg_partons.begin();
      if (const auto psgen_partons = partons();  // ensure parton system is compatible
          static_cast<int>(mg_hadron_parton) != static_cast<int>(psgen_partons.at(electron_first_ ? 1 : 0)))
        throw CG_FATAL("ElectronPartonProcessBuilder")
            << "MadGraph unpacked process incoming state (" << process().intermediatePartons()
            << ") is incompatible with user-steered incoming fluxes particles (" << psgen_partons << ").";
      prepareSteeringCard();
    }

  private:
    Momentum leptonBeam(bool positive) const {
      const auto e = 0.5 * std::sqrt(s());
      return Momentum::fromPxPyPzE(0., 0., positive ? e : -e, e);
    }

    double matrixElementNormalisation() const override {
      const bool lep_pos = electron_first_;
      const auto s_eg =
          (lep_pos ? leptonBeam(true) + q2() : leptonBeam(false) + q1()).mass2();  // beam lepton + hadron-side photon
      const auto x_had = lep_pos ? x2() : x1(), x_lep = lep_pos ? x1() : x2();
      if (s_eg <= 0. || x_had <= 0. || x_lep >= 1.)
        return 0.;
      const auto lep_q2 = std::fabs((lep_pos ? leptonBeam(true) - pX() : leptonBeam(false) - pY()).mass2());
      if (!kinematics()
               .cuts()
               .initial.q2
               .at(static_cast<size_t>(electron_first_ ? PartonSampler::Side::positive : PartonSampler::Side::negative))
               .contains(lep_q2))  // lepton-side photon virtuality not sampled: apply parton virtuality cut explicitly
        return 0.;
      const auto norm = 1. / (2. * std::pow(2. * M_PI, 3) * x_had * s() * s_eg * (1. - x_lep));
      CG_DEBUG_LOOP("ElectronPartonProcessBuilder:norm")
          << "x_had=" << x_had << ", x_lep=" << x_lep << ", s_eg=" << s_eg << ", lepton q2=" << lep_q2
          << ", norm=" << norm << ".";
      return norm;
    }

    void fillMomenta() override {
      Momentum central_sum;
      for (size_t j = 0, n = central().size(); j < n; ++j)
        central_sum += pc(j);
      if (electron_first_)
        pX() = leptonBeam(true) + q2() - central_sum;
      else
        pY() = leptonBeam(false) + q1() - central_sum;
      size_t i = 0;
      if (electron_first_) {
        process().setMomentum(i++, leptonBeam(true));  // first "parton": beam electron
        process().setMomentum(i++, q2());              // second "parton": parton-from-hadron
        process().setMomentum(i++, pX());
        for (size_t j = 0, n = central().size(); j < n; ++j)
          process().setMomentum(i++, pc(j));  // outgoing central particles
      } else {
        process().setMomentum(i++, q1());               // first "parton": parton-from-hadron
        process().setMomentum(i++, leptonBeam(false));  // second "parton": beam electron
        for (size_t j = 0, n = central().size(); j < n; ++j)
          process().setMomentum(i++, pc(j));  // outgoing central particles
        process().setMomentum(i, pY());
      }
    }

    bool electron_first_{true};  ///< is the lepton the first (or last) particle of the MadGraph central system?
  };
}  // namespace cepgen::mg5amc
using MadGraphElectronPartonProcessBuilder = cepgen::mg5amc::ElectronPartonProcessBuilder;
REGISTER_PROCESS("mg5_aMC:eh", MadGraphElectronPartonProcessBuilder);
