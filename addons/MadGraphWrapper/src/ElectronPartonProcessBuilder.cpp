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

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGenMadGraph/Interface.h"
#include "CepGenMadGraph/Process.h"
#include "CepGenMadGraph/ProcessFactory.h"
#include "CepGenMadGraph/Utils.h"

using namespace cepgen;
using namespace std::string_literals;

namespace cepgen::mg5amc {
  class ElectronPartonProcessBuilder : public proc::FactorisedProcess {
  public:
    explicit ElectronPartonProcessBuilder(const ParametersList& params, bool load_library = true)
        : FactorisedProcess(params, {}) {
      if (load_library)
        loadMG5Library();
      CG_DEBUG("ElectronPartonProcessBuilder")
          << "List of MadGraph process registered in the runtime database: " << ProcessFactory::get().modules() << ".";
      // once MadGraph process library is loaded into runtime environment, can define its wrapper object
      mg5_proc_ = ProcessFactory::get().build(normalise(steer<std::string>("process")));
      if (mg5_proc_->centralSystem().empty())
        throw CG_FATAL("ElectronPartonProcessBuilder")
            << "Failed to retrieve produced particles system from MadGraph process:\n"
            << mg5_proc_->description().validate(mg5_proc_->parameters()) << ".";
      auto central_system = mg5_proc_->centralSystem();
      if (std::abs(*central_system.begin()) == PDG::electron) {
        e_pdg_ = *central_system.begin();
        mode_ = Mode::electron_parton;
        central_system.erase(central_system.begin());
      } else if (std::abs(*central_system.rbegin()) == PDG::electron) {
        e_pdg_ = *central_system.rbegin();
        mode_ = Mode::parton_electron;
        central_system.pop_back();
      } else
        throw CG_FATAL("ElectronPartonProcessBuilder")
            << "No electron/positron found in mg5_aMC process particles list.";
      phase_space_generator_->setCentral(central_system);
    }

    proc::ProcessPtr clone() const override {
      return std::make_unique<ElectronPartonProcessBuilder>(parameters(), false);
    }

    void addEventContent() override {
      const auto central_system = phase_space_generator_->central();
      Process::setEventContent(
          {{Particle::Role::IncomingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
           {Particle::Role::IncomingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
           {Particle::Role::OutgoingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
           {Particle::Role::OutgoingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
           {Particle::Role::CentralSystem, spdgids_t(central_system.begin(), central_system.end())}});
    }

    static ParametersDescription description() {
      auto desc = FactorisedProcess::description();
      desc.setDescription("MadGraph_aMC electron-parton process builder");
      desc.add("lib", ""s).setDescription("Precompiled library for this process definition");
      desc.add("parametersCard", "param_card.dat"s).setDescription("Runtime MadGraph parameters card");
      desc += Interface::description();
      return desc;
    }

    void prepareFactorisedPhaseSpace() override {
      if (const auto psgen_partons = phase_space_generator_->partons();
          (mode_ == Mode::parton_electron &&
           *mg5_proc_->intermediatePartons().begin() != static_cast<int>(psgen_partons.at(0))) ||
          (mode_ == Mode::electron_parton &&
           *mg5_proc_->intermediatePartons().rbegin() != static_cast<int>(psgen_partons.at(1))))
        throw CG_FATAL("ElectronPartonProcessBuilder")
            << "MadGraph unpacked process incoming state (" << mg5_proc_->intermediatePartons()
            << ") is incompatible with user-steered incoming fluxes particles (" << phase_space_generator_->partons()
            << ").";
      if (const auto params_card = steer<std::string>("parametersCard"); !params_card.empty()) {
        CG_INFO("ElectronPartonProcessBuilder")
            << "Preparing process kinematics for card at \"" << params_card << "\".";
        const auto unsteered_pcard = Interface::extractParamCardParameters(utils::readFile(params_card));
        CG_DEBUG("ElectronPartonProcessBuilder") << "Unsteered parameters card:\n" << unsteered_pcard;
        if (const auto mod_params = steer<ParametersList>("modelParameters"); !mod_params.empty()) {
          const auto steered_pcard = unsteered_pcard.steer(mod_params);
          CG_DEBUG("ElectronPartonProcessBuilder") << "User-steered parameters:" << mod_params << "\n"
                                                   << "Steered parameters card:\n"
                                                   << steered_pcard;
          std::ofstream params_card_steered(params_card);
          params_card_steered << Interface::generateParamCard(steered_pcard);
          params_card_steered.close();
        }
        mg5_proc_->initialise(params_card);
      }
    }
    double computeFactorisedMatrixElement() override {
      if (!mg5_proc_)
        throw CG_FATAL("ElectronPartonProcessBuilder:eval") << "Process not properly linked!";
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
        mg5_proc_->setMomentum(i++, pA());  // first "parton" is beam electron
        mg5_proc_->setMomentum(i++, q2());  // second "parton" is parton-from-hadron
        mg5_proc_->setMomentum(i++, pX());
        for (size_t j = 0; j < phase_space_generator_->central().size(); ++j)
          mg5_proc_->setMomentum(i++, pc(j));  // outgoing central particles
      } else if (mode_ == Mode::parton_electron) {
        size_t i = 0;
        mg5_proc_->setMomentum(i++, q1());  // first "parton" is parton-from-hadron
        mg5_proc_->setMomentum(i++, pB());  // second "parton" is beam electron
        for (size_t j = 0; j < phase_space_generator_->central().size(); ++j)
          mg5_proc_->setMomentum(i++, pc(j));  // outgoing central particles
        mg5_proc_->setMomentum(i++, pY());
      } else
        throw CG_FATAL("ElectronPartonProcessBuilder:eval") << "Invalid beams mode: " << mode_ << ".";
      if (const auto weight = mg5_proc_->eval(); utils::positive(weight))
        return weight * std::pow(shat(), -2);
      return 0.;
    }

  private:
    void loadMG5Library() const {
      utils::AbortHandler();
      try {
        if (const auto& lib_file = steer<std::string>("lib"); !lib_file.empty())  // user-provided library file
          loadLibrary(lib_file);
        else {  // library has to be generated from mg5_aMC directives
          const Interface interface(params_);
          loadLibrary(interface.run());
        }
      } catch (const utils::RunAbortedException&) {
        CG_FATAL("ElectronPartonProcessBuilder") << "MadGraph_aMC process generation aborted.";
      }
    }
    std::unique_ptr<mg5amc::Process> mg5_proc_;
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
