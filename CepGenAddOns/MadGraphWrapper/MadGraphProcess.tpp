/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2023  Laurent Forthomme
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

//=============================================================================
// NOLI SE TANGERE
#include "CPPProcess.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Math.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

using namespace cepgen;

MadGraphProcess::MadGraphProcess()
    : proc_(new CPPProcess),
      name_("XXX_PROC_NAME_XXX"),
      descr_("XXX_PROC_DESCRIPTION_XXX"),
      incoming_pdgids_({XXX_PART1_XXX, XXX_PART2_XXX}),
      central_pdgids_({XXX_OUT_PART_XXX}) {
  CG_INFO("MadGraphProcess") << "Process considered: " << proc_->name() << ". "
                             << "Incoming particles: " << incoming_pdgids_ << ", outgoing system: " << central_pdgids_
                             << ".";
}

MadGraphProcess::~MadGraphProcess() = default;

void MadGraphProcess::initialise(const std::string& param_card) {
  try {
    proc_->initProc(param_card);
  } catch (const char* chr) {
    throw CG_FATAL("MadGraphProcess:init")
        << "Failed to initialise parameters card at \"" << param_card << "\":\n\t" << chr;
  }
  if (proc_->nprocesses > 1)
    throw CG_FATAL("MadGraphProcess:init") << "Multi-processes matrix elements are not (yet) supported!";
  if (proc_->ninitial != 2)
    throw CG_FATAL("MadGraphProcess:init") << "Currently only 2->N processes are supported!";

  CG_DEBUG("MadGraphProcess:init") << "External particles masses (partons + central system): " << proc_->getMasses()
                                   << ".";

  mom_.clear();
  for (size_t i = 0; i < proc_->nexternal; ++i)
    mom_.emplace_back(new double[4]{proc_->getMasses().at(i), 0., 0., 0.});
  momenta_.reserve(proc_->nexternal);
}

double MadGraphProcess::eval() {
  proc_->setMomenta(mom_);
  proc_->sigmaKin();
  const double* me = proc_->getMatrixElements();
  if (!utils::positive(me[0]))
    return 0.;

  CG_DEBUG_LOOP("MadGraphProcess:eval").log([&](auto& log) {
    log << "Dump of event kinematics\n\t"
        << "Incoming partons 4-momenta:   " << std::vector<double>(mom_[0], mom_[0] + 4) << ", "
        << std::vector<double>(mom_[1], mom_[1] + 4) << "\n\t"
        << "Outgoing particles 4-momenta: ";
    std::string sep;
    for (size_t i = 0; i < proc_->nexternal - 2; ++i)
      log << sep << std::vector<double>(mom_[i + 2], mom_[i + 2] + 4), sep = ", ";
    log << "\n\tResulting matrix element: " << me[0] << ".";
  });
  return me[0];
}

const std::vector<Momentum>& MadGraphProcess::momenta() {
  const auto& p4 = proc_->getMomenta();
  // cast it to the member attribute and return it
  for (size_t i = 0; i < p4.size(); ++i)
    momenta_[i] = Momentum::fromPxPyPzE(p4[i][1], p4[i][2], p4[i][3], p4[i][0]);
  return momenta_;
}
//=============================================================================
