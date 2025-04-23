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
  CG_DEBUG("ProcessBuilder") << "List of MadGraph process registered in the runtime database: "
                             << ProcessFactory::get().modules() << ".";
  // once MadGraph process library is loaded into runtime environment, can define its wrapper object
  mg5_proc_ = ProcessFactory::get().build(normalise(steer<std::string>("process")));
  if (mg5_proc_->centralSystem().empty())
    throw CG_FATAL("ProcessBuilder") << "Failed to retrieve produced particles system from MadGraph process:\n"
                                     << mg5_proc_->description().validate(mg5_proc_->parameters()) << ".";
}

void ProcessBuilder::addEventContent() {
  const auto central_system = phase_space_generator_->central();
  setEventContent({{Particle::Role::IncomingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                   {Particle::Role::IncomingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                   {Particle::Role::OutgoingBeam1, {kinematics().incomingBeams().positive().integerPdgId()}},
                   {Particle::Role::OutgoingBeam2, {kinematics().incomingBeams().negative().integerPdgId()}},
                   {Particle::Role::CentralSystem, spdgids_t(central_system.begin(), central_system.end())}});
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
    CG_FATAL("ProcessBuilder") << "MadGraph_aMC process generation aborted.";
  }
}
