/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2024  Laurent Forthomme
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
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcessFactory.h"
#include "CepGenAddOns/MadGraphWrapper/Utils.h"

using namespace cepgen;

class MadGraphProcessBuilder : public proc::FactorisedProcess {
public:
  explicit MadGraphProcessBuilder(const ParametersList& params, bool load_library = true)
      : FactorisedProcess(params, {}) {
    if (load_library)
      loadMG5Library();
    CG_DEBUG("MadGraphProcessBuilder") << "List of MadGraph process registered in the runtime database: "
                                       << MadGraphProcessFactory::get().modules() << ".";
    // once MadGraph process library is loaded into runtime environment, can define its wrapper object
    mg5_proc_ = MadGraphProcessFactory::get().build(mg5amc::normalise(steer<std::string>("process")));
    if (mg5_proc_->centralSystem().empty())
      throw CG_FATAL("MadGraphProcessBuilder")
          << "Failed to retrieve produced particles system from MadGraph process:\n"
          << mg5_proc_->description().validate(mg5_proc_->parameters()) << ".";
    psgen_->setCentral(mg5_proc_->centralSystem());
  }

  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new MadGraphProcessBuilder(parameters(), false)); }

  void addEventContent() override {
    const auto mg5_proc_cent = mg5_proc_->centralSystem();
    Process::setEventContent({{Particle::IncomingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                              {Particle::IncomingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                              {Particle::OutgoingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                              {Particle::OutgoingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                              {Particle::CentralSystem, spdgids_t(mg5_proc_cent.begin(), mg5_proc_cent.end())}});
  }

  static ParametersDescription description() {
    auto desc = FactorisedProcess::description();
    desc.setDescription("MadGraph_aMC process builder");
    desc.add<std::string>("lib", "").setDescription("Precompiled library for this process definition");
    desc.add<std::string>("parametersCard", "param_card.dat").setDescription("Runtime MadGraph parameters card");
    desc += MadGraphInterface::description();
    return desc;
  }

  void prepareFactorisedPhaseSpace() override {
    const auto psgen_partons = psgen_->partons();
    if (mg5_proc_->intermediatePartons() != std::vector<int>(psgen_partons.begin(), psgen_partons.end()))
      throw CG_FATAL("MadGraphProcessBuilder")
          << "MadGraph unpacked process incoming state (" << mg5_proc_->intermediatePartons() << ") "
          << "is incompatible with user-steered incoming fluxes particles (" << psgen_->partons() << ").";
    if (const auto params_card = steer<std::string>("parametersCard"); !params_card.empty()) {
      CG_INFO("MadGraphProcessBuilder") << "Preparing process kinematics for card at \"" << params_card << "\".";
      if (const auto mod_params = steer<ParametersList>("modelParameters"); !mod_params.empty()) {
        const auto unsteered_pcard_txt = utils::readFile(params_card);
        const auto steered_pcard = MadGraphInterface::extractParamCardParameters(unsteered_pcard_txt).steer(mod_params);
        CG_DEBUG("MadGraphProcessBuilder") << "Unsteered parameters card:\n"
                                           << unsteered_pcard_txt << "\n\n"
                                           << std::string(50, '-') << "\n"
                                           << "Steered parameters card:\n"
                                           << steered_pcard;
        std::ofstream params_card_steered(params_card);
        params_card_steered << MadGraphInterface::generateParamCard(steered_pcard);
        params_card_steered.close();
      }
      mg5_proc_->initialise(params_card);
    }
  }
  double computeFactorisedMatrixElement() override {
    if (!mg5_proc_)
      CG_FATAL("MadGraphProcessBuilder:eval") << "Process not properly linked!";
    if (!kinematics().cuts().initial.contain(event()(Particle::Role::Parton1)) ||
        !kinematics().cuts().initial.contain(event()(Particle::Role::Parton2)))
      return 0.;
    if (!kinematics().cuts().central.contain(event()(Particle::Role::CentralSystem)))
      return 0.;

    CG_DEBUG_LOOP("MadGraphProcessBuilder:eval")
        << "Particles content:\n"
        << "incoming: " << q1() << " (m=" << q1().mass() << "), " << q2() << " (m=" << q2().mass() << ")\n"
        << "outgoing: " << pc(0) << " (m=" << pc(0).mass() << "), " << pc(1) << " (m=" << pc(1).mass() << ").";
    mg5_proc_->setMomentum(0, q1());   // first incoming parton
    mg5_proc_->setMomentum(1, q2());   // second incoming parton
    mg5_proc_->setMomentum(2, pc(0));  // first outgoing central particle
    mg5_proc_->setMomentum(3, pc(1));  // second outgoing central particle
    if (const auto weight = mg5_proc_->eval(); utils::positive(weight))
      return weight * std::pow(shat(), -2);
    return 0.;
  }

private:
  void loadMG5Library() {
    utils::AbortHandler();
    try {
      if (const auto& lib_file = steer<std::string>("lib"); !lib_file.empty())
        loadLibrary(lib_file);
      else {
        const MadGraphInterface interf(params_);
        loadLibrary(interf.run());
      }
    } catch (const utils::RunAbortedException&) {
      CG_FATAL("MadGraphProcessBuilder") << "MadGraph_aMC process generation aborted.";
    }
  }

  std::unique_ptr<MadGraphProcess> mg5_proc_;
};
REGISTER_PROCESS("mg5_aMC", MadGraphProcessBuilder);
