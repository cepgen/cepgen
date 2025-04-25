/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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
        mode_ = Mode::electron_parton;
        e_pdg_ = *central_system.begin();              // electron, or positron?
        central_system.erase(central_system.begin());  // first particle is lepton
      } else if (std::abs(*central_system.rbegin()) == PDG::electron) {
        mode_ = Mode::parton_electron;
        e_pdg_ = *central_system.rbegin();  // electron, or positron?
        central_system.pop_back();          // last particle is lepton
      } else
        throw CG_FATAL("ElectronPartonProcessBuilder")
            << "No electron/positron found in mg5_aMC process particles list.";
      phase_space_generator_->setCentral(central_system);  // electron/positron stripped off central system
    }

    proc::ProcessPtr clone() const override {
      return std::make_unique<ElectronPartonProcessBuilder>(parameters(), false);
    }

    static ParametersDescription description() {
      auto desc = ProcessBuilder::description();
      desc.setDescription("MadGraph_aMC electron-parton process builder");
      return desc;
    }

    void prepareFactorisedPhaseSpace() override {
      if (const auto psgen_partons = phase_space_generator_->partons();  // ensure parton system is compatible
          (mode_ == Mode::parton_electron &&
           *process().intermediatePartons().begin() != static_cast<int>(psgen_partons.at(0))) ||
          (mode_ == Mode::electron_parton &&
           *process().intermediatePartons().rbegin() != static_cast<int>(psgen_partons.at(1))))
        throw CG_FATAL("ElectronPartonProcessBuilder")
            << "MadGraph unpacked process incoming state (" << process().intermediatePartons()
            << ") is incompatible with user-steered incoming fluxes particles (" << phase_space_generator_->partons()
            << ").";
      prepareSteeringCard();
    }
    double computeFactorisedMatrixElement() override {
      if (!kinematics().cuts().initial.contain(event()(Particle::Role::Parton1)) ||
          !kinematics().cuts().initial.contain(event()(Particle::Role::Parton2)))
        return 0.;
      if (!kinematics().cuts().central.contain(event()(Particle::Role::CentralSystem)))
        return 0.;
      pX() = (pA() - q1()).setMass2(mX2());
      pY() = (pB() - q2()).setMass2(mY2());
      CG_DEBUG_LOOP("ElectronPartonProcessBuilder:eval")
          << "Particles content:\n"
          << "incoming: " << q1() << " (m=" << q1().mass() << "), " << q2() << " (m=" << q2().mass() << ")\n"
          << "outgoing: " << pc(0) << " (m=" << pc(0).mass() << "), " << pc(1) << " (m=" << pc(1).mass() << ").";
      if (mode_ == Mode::electron_parton) {
        size_t i = 0;
        process().setMomentum(i++, pA());  // first "parton": beam electron
        process().setMomentum(i++, q2());  // second "parton": parton-from-hadron
        process().setMomentum(i++, pX());
        for (size_t j = 0; j < phase_space_generator_->central().size(); ++j)
          process().setMomentum(i++, pc(j));  // outgoing central particles
      } else if (mode_ == Mode::parton_electron) {
        size_t i = 0;
        process().setMomentum(i++, q1());  // first "parton": parton-from-hadron
        process().setMomentum(i++, pB());  // second "parton": beam electron
        for (size_t j = 0; j < phase_space_generator_->central().size(); ++j)
          process().setMomentum(i++, pc(j));  // outgoing central particles
        process().setMomentum(i, pY());
      } else
        throw CG_FATAL("ElectronPartonProcessBuilder:eval") << "Invalid beams mode: " << mode_ << ".";
      if (const auto weight = process().eval(); utils::positive(weight))
        return weight * std::pow(shat(), -2);
      return 0.;
    }

  private:
    spdgid_t e_pdg_{11};
    enum struct Mode { electron_parton, parton_electron } mode_{Mode::electron_parton};
    friend std::ostream& operator<<(std::ostream& os, const Mode& mode) {
      switch (mode) {
        case Mode::electron_parton:
          return os << "electron-parton";
        case Mode::parton_electron:
          return os << "parton-electron";
      }
      return os << "invalid";
    }
  };
}  // namespace cepgen::mg5amc
using MadGraphElectronPartonProcessBuilder = cepgen::mg5amc::ElectronPartonProcessBuilder;
REGISTER_PROCESS("mg5_aMC:eh", MadGraphElectronPartonProcessBuilder);
