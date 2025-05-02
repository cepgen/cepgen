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
#include "CepGen/Generator.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGenMadGraph/Interface.h"
#include "CepGenMadGraph/Process.h"
#include "CepGenMadGraph/ProcessBuilder.h"
#include "CepGenMadGraph/ProcessFactory.h"
#include "CepGenMadGraph/Utils.h"

using namespace cepgen;
using namespace cepgen::mg5amc;
using namespace std::string_literals;

ProcessBuilder::ProcessBuilder(const ParametersList& params, bool load_library) : FactorisedProcess(params, {}) {
  if (load_library)
    loadMG5Library();
  CG_DEBUG("mg5amc:ProcessBuilder") << "List of MadGraph process registered in the runtime database: "
                                    << ProcessFactory::get().modules() << ".";
  // once MadGraph process library is loaded into runtime environment, can define its wrapper object
  mg5_proc_ = ProcessFactory::get().build(normalise(steer<std::string>("process")));
  if (mg5_proc_->centralSystem().empty())
    throw CG_FATAL("mg5amc:ProcessBuilder") << "Failed to retrieve produced particles system from MadGraph process:\n"
                                            << mg5_proc_->description().validate(mg5_proc_->parameters()) << ".";
}

ParametersDescription ProcessBuilder::description() {
  auto desc = FactorisedProcess::description();
  desc.setDescription("MadGraph_aMC process builder");
  desc.add("lib", ""s).setDescription("Precompiled library for this process definition");
  desc.add("parametersCard", "param_card.dat"s).setDescription("Runtime MadGraph parameters card");
  desc += Interface::description();
  return desc;
}

void ProcessBuilder::loadMG5Library() const {
  utils::AbortHandler();
  try {
    if (const auto& lib_file = steer<std::string>("lib"); !lib_file.empty())  // user-provided library file
      loadLibrary(lib_file);
    else {  // library has to be generated from mg5_aMC directives
      const Interface interface(params_);
      loadLibrary(interface.run());
    }
  } catch (const utils::RunAbortedException&) {
    CG_FATAL("mg5amc:ProcessBuilder") << "MadGraph_aMC process generation aborted.";
  }
}

void ProcessBuilder::prepareSteeringCard() const {
  if (const auto params_card = steer<std::string>("parametersCard"); !params_card.empty()) {
    CG_INFO("mg5amc:ProcessBuilder") << "Preparing process kinematics for card at \"" << params_card << "\".";
    const auto unsteered_pcard = Interface::extractParamCardParameters(utils::readFile(params_card));
    CG_DEBUG("mg5amc:ProcessBuilder") << "Unsteered parameters card:\n" << unsteered_pcard;
    if (const auto mod_params = steer<ParametersList>("modelParameters"); !mod_params.empty()) {
      const auto steered_pcard = unsteered_pcard.steer(mod_params);
      CG_DEBUG("mg5amc:ProcessBuilder") << "User-steered parameters:" << mod_params << "\n"
                                        << "Steered parameters card:\n"
                                        << steered_pcard;
      std::ofstream params_card_steered(params_card);
      params_card_steered << Interface::generateParamCard(steered_pcard);
      params_card_steered.close();
    }
    mg5_proc_->initialise(params_card);
  }
}

mg5amc::Process& ProcessBuilder::process() const {
  if (!mg5_proc_)
    CG_FATAL("mg5amc:GeneralProcessBuilder:eval") << "Process not properly linked!";
  return *mg5_proc_;
}
