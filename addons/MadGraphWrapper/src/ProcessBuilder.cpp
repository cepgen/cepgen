/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2026  Laurent Forthomme
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

#include <CepGen/Core/Exception.h>
#include <CepGen/Event/Event.h>
#include <CepGen/Generator.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Process/FactorisedProcess.h>
#include <CepGen/Utils/AbortHandler.h>
#include <CepGen/Utils/Math.h>
#include <CepGen/Utils/RandomGenerator.h>

#include <fstream>

#include "CepGenMadGraph/Interface.h"
#include "CepGenMadGraph/Process.h"
#include "CepGenMadGraph/ProcessBuilder.h"
#include "CepGenMadGraph/ProcessFactory.h"
#include "CepGenMadGraph/Utils.h"

using namespace cepgen;
using namespace cepgen::mg5amc;
using namespace std::string_literals;

ProcessBuilder::ProcessBuilder(const ParametersList& params, bool load_library)
    : FactorisedProcess(params), owns_library_(load_library), library_filename_(steer<std::string>("lib")) {
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

ProcessBuilder::~ProcessBuilder() {
  mg5_proc_.reset();  // call the destructor, as the library will be invalidated
  if (!owns_library_)
    return;  // clones share the library of the process they were built from
  if (!library_filename_.empty() && !unloadLibrary(library_filename_))
    CG_ERROR("mg5amc:~ProcessBuilder") << "Failed to unload library '" << library_filename_ << "'.";
  if (steer<bool>("removeLibrary")) {
    std::error_code error_code;
    if (library_filename_.empty() || !std::filesystem::remove_all(library_filename_, error_code))
      CG_WARNING("mg5amc:~ProcessBuilder")
          << "Failed to remove generated library '" << library_filename_ << "'. Error: " << error_code.message() << ".";
    else
      CG_DEBUG("mg5amc:~ProcessBuilder") << "Successfully removed generated library '" << library_filename_ << "'.";
  }
}

ParametersDescription ProcessBuilder::description() {
  auto desc = FactorisedProcess::description();
  desc.setDescription("MadGraph_aMC process builder");
  desc.add("lib", ""s).setDescription("Precompiled library for this process definition");
  desc.add("parametersCard", "param_card.dat"s).setDescription("Runtime MadGraph parameters card");
  desc.add("removeLibrary", false).setDescription("Remove the library after usage?");
  desc += Interface::description();
  return desc;
}

void ProcessBuilder::loadMG5Library() {
  utils::AbortHandler();
  try {
    if (library_filename_.empty()) {  // library has to be generated from mg5_aMC directives
      const Interface interface(params_);
      library_filename_ = interface.run();
      // record the process location in the steering parameters, so that it is inherited by all process clones
      params_.set("tmpDir", fs::path(library_filename_).parent_path().string());
    }
    if (!loadLibrary(library_filename_))
      throw CG_ERROR("mg5amc:ProcessBuilder:loadMG5Library")
          << "Failed to load MadGraph_aMC process library '" << library_filename_ << "'.";
    CG_DEBUG("mg5amc:ProcessBuilder:loadMG5Library")
        << "MadGraph_aMC process loaded from library '" << library_filename_ << "'.";
  } catch (const utils::RunAbortedException&) {
    throw CG_FATAL("mg5amc:ProcessBuilder:loadMG5Library") << "MadGraph_aMC process generation aborted.";
  }
}

void ProcessBuilder::prepareSteeringCard() const {
  if (const auto params_card = parametersCardPath(); !params_card.empty()) {
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

double ProcessBuilder::computeFactorisedMatrixElement() {
  if (!kinematics().cuts().initial.contain(event()(Particle::Role::Parton1)) ||
      !kinematics().cuts().initial.contain(event()(Particle::Role::Parton2)) ||
      !kinematics().cuts().central.contain(event()(Particle::Role::CentralSystem)))
    return 0.;
  CG_DEBUG_LOOP("mg5amc:ProcessBuilder:eval")
      << "Incoming partons: " << q1() << " (m=" << q1().mass() << "), " << q2() << " (m=" << q2().mass() << ").";
  fillMomenta();
  if (const auto weight = process().eval(); utils::positive(weight)) {
    const auto s = shat();
    return weight / (s * s);
  }
  return 0.;
}

std::string ProcessBuilder::parametersCardPath() const {
  const fs::path card = steer<std::string>("parametersCard");
  if (card.empty() || (fs::exists(card) && !fs::is_symlink(card)))
    return card.string();  // trivial case: cards were found in lib directory or with user-provided path
  const auto& lib_path = steer<std::string>("lib");
  const auto& standalone_path = steer<std::string>("standaloneCppPath");
  const fs::path process_path = !lib_path.empty()          ? fs::path{lib_path}.parent_path()
                                : !standalone_path.empty() ? fs::path{standalone_path}
                                                           : Interface::temporaryPath(steer<std::string>("tmpDir"));
  if (const auto generated_card = process_path / "Cards" / card; fs::exists(generated_card))
    return generated_card.string();
  if (fs::exists(card))
    return card.string();
  throw CG_FATAL("mg5amc:ProcessBuilder") << "Failed to find the parameters card '" << card
                                          << "' in the working directory, or in " << process_path / "Cards" << ".";
}
