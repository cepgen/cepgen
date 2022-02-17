/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Processes/Process2to4.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

using namespace cepgen;

class MadGraphProcessBuilder : public proc::Process2to4 {
public:
  MadGraphProcessBuilder(const ParametersList&);
  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new MadGraphProcessBuilder(*this)); }

  static ParametersDescription description();

  void prepareProcessKinematics() override;
  double computeCentralMatrixElement() const override;

private:
  std::shared_ptr<MadGraphProcess> mg5_proc_;
};

MadGraphProcessBuilder::MadGraphProcessBuilder(const ParametersList& params)
    : Process2to4(params, std::array<pdgid_t, 2>{}, 0) {
  utils::AbortHandler();
  try {
    if (params_.has<std::string>("lib"))
      loadLibrary(steer<std::string>("lib"));
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
  setIntermediatePartons({(pdgid_t)interm_part[0], (pdgid_t)interm_part[1]});
  setProducedParticles(std::vector<pdgid_t>(cent_sys.begin(), cent_sys.end()));
}

void MadGraphProcessBuilder::prepareProcessKinematics() {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  mg5_proc_->initialise(steer<std::string>("parametersCard"));
}

double MadGraphProcessBuilder::computeCentralMatrixElement() const {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  mg5_proc_->setMomentum(0, q1_);    // first incoming parton
  mg5_proc_->setMomentum(1, q2_);    // second incoming parton
  mg5_proc_->setMomentum(2, p_c1_);  // first outgoing central particle
  mg5_proc_->setMomentum(3, p_c2_);  // second outgoing central particle

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

REGISTER_PROCESS("mg5_aMC", MadGraphProcessBuilder)
