/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2022  Laurent Forthomme
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
#include "CepGen/Process/Process2to4.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

using namespace cepgen;

class MadGraphProcessBuilder : public proc::Process2to4 {
public:
  explicit MadGraphProcessBuilder(const ParametersList&);
  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new MadGraphProcessBuilder(*this)); }

  static ParametersDescription description();

  void prepareProcessKinematics() override;
  double computeCentralMatrixElement() const override;

private:
  std::shared_ptr<MadGraphProcess> mg5_proc_;
};

MadGraphProcessBuilder::MadGraphProcessBuilder(const ParametersList& params) : Process2to4(params, {}) {
  utils::AbortHandler();
  try {
    const auto& lib_file = steer<std::string>("lib");
    if (!lib_file.empty())
      loadLibrary(lib_file);
    else {
      const MadGraphInterface interf(params);
      loadLibrary(interf.run());
    }
  } catch (const utils::RunAbortedException&) {
    CG_FATAL("MadGraphProcessBuilder") << "MadGraph_aMC process generation aborted.";
  }
  // once MadGraph process library is loaded into runtime environment, can define its wrapper object
  mg5_proc_.reset(new MadGraphProcess);
  const auto& interm_part = mg5_proc_->intermediatePartons();
  const auto& cent_sys = mg5_proc_->centralSystem();
  CG_DEBUG("MadGraphProcessBuilder") << "MadGraph_aMC process created for:\n\t"
                                     << "* interm. parts.: " << interm_part << "\n\t"
                                     << "* central system: " << cent_sys << ".";
  kinematics().incomingBeams().positive().setPdgId(interm_part[0]);
  kinematics().incomingBeams().negative().setPdgId(interm_part[1]);
  produced_parts_ = std::vector<pdgid_t>(cent_sys.begin(), cent_sys.end());
}

void MadGraphProcessBuilder::prepareProcessKinematics() {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  const auto& params_card = steer<std::string>("parametersCard");
  CG_INFO("MadGRaphProcessBuilder") << "Preparing process kinematics for card at \"" << params_card << "\".";
  mg5_proc_->initialise(params_card);
}

double MadGraphProcessBuilder::computeCentralMatrixElement() const {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder:eval") << "Process not properly linked!";

  CG_DEBUG_LOOP("MadGraphProcessBuilder:eval")
      << "Particles content:\n"
      << "incoming: " << q1() << " (m=" << q1().mass() << "), " << q2() << " (m=" << q2().mass() << ")\n"
      << "outgoing: " << pc(0) << " (m=" << pc(0).mass() << "), " << pc(1) << " (m=" << pc(1).mass() << ").";
  mg5_proc_->setMomentum(0, q1());   // first incoming parton
  mg5_proc_->setMomentum(1, q2());   // second incoming parton
  mg5_proc_->setMomentum(2, pc(0));  // first outgoing central particle
  mg5_proc_->setMomentum(3, pc(1));  // second outgoing central particle

  return mg5_proc_->eval();
}

ParametersDescription MadGraphProcessBuilder::description() {
  auto desc = Process2to4::description();
  desc.setDescription("MadGraph_aMC process builder");
  desc.add<std::string>("lib", "").setDescription("Precompiled library for this process definition");
  desc.add<std::string>("parametersCard", "param_card.dat").setDescription("Runtime MadGraph parameters card");
  desc += MadGraphInterface::description();
  return desc;
}

REGISTER_PROCESS("mg5_aMC", MadGraphProcessBuilder);
