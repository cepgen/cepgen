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
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

using namespace cepgen;

class MadGraphProcessBuilder : public proc::Process2to4 {
public:
  MadGraphProcessBuilder(const ParametersList& params);
  proc::ProcessPtr clone() const override { return proc::ProcessPtr(new MadGraphProcessBuilder(*this)); }
  static std::string description() { return "MadGraph_aMC process builder"; }

  void prepareProcessKinematics() override;
  double computeCentralMatrixElement() const override;

  /*void preparePhaseSpace() override;
    double computeKTFactorisedMatrixElement() override;
    void fillCentralParticlesKinematics() override;*/

private:
  std::shared_ptr<MadGraphProcess> mg5_proc_;
};

extern std::string madgraph_process_name();

MadGraphProcessBuilder::MadGraphProcessBuilder(const ParametersList& params)
    :  //KTProcess( params, std::array<pdgid_t,2>{}, std::vector<pdgid_t>{} )
      Process2to4(params, std::array<pdgid_t, 2>{}, 0) {
  if (params.has<std::string>("lib"))
    loadLibrary(params.get<std::string>("lib"));
  else {
    const MadGraphInterface interf(params);
    loadLibrary(interf.run());
  }
  //--- once MadGraph process library is loaded into runtime environment
  //    can define its wrapper object
  mg5_proc_.reset(new MadGraphProcess);
}

void
//MadGraphProcessBuilder::preparePhaseSpace()
MadGraphProcessBuilder::prepareProcessKinematics() {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  mg5_proc_->initialise(params_.get<std::string>("parametersCard", "param_card.dat"));
}

/*void MadGraphProcessBuilder::fillCentralParticlesKinematics() {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  const auto& parts = mg5_proc_->momenta();
}*/

double
//MadGraphProcessBuilder::computeKTFactorisedMatrixElement()
MadGraphProcessBuilder::computeCentralMatrixElement() const {
  if (!mg5_proc_)
    CG_FATAL("MadGraphProcessBuilder") << "Process not properly linked!";

  mg5_proc_->setMomentum(0, q1_);  // first incoming parton
  mg5_proc_->setMomentum(1, q2_);  // second incoming parton
  mg5_proc_->setMomentum(2, p_c1_);
  mg5_proc_->setMomentum(3, p_c2_);

  return mg5_proc_->eval();
}

REGISTER_PROCESS("mg5_aMC", MadGraphProcessBuilder)
