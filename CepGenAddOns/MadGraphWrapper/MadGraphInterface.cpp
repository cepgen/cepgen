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

#ifndef MADGRAPH_BIN
#error "*** MADGRAPH_BIN variable not set! ***"
#endif
#ifndef MADGRAPH_PROC_TMPL
#error "*** MADGRAPH_PROC_TMPL variable not set! ***"
#endif
#ifndef CC_CFLAGS
#error "*** CC_CFLAGS variable not set! ***"
#endif

#include <array>
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Caller.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"
#include "CepGenAddOns/MadGraphWrapper/Utils.h"
#include "CepGenAddOns/PythonWrapper/Environment.h"

namespace cepgen {
  std::unordered_map<std::string, pdgid_t> MadGraphInterface::mg5_parts_ = {
      {"d", (pdgid_t)1},     {"d~", (pdgid_t)1},    {"u", (pdgid_t)2},    {"u~", (pdgid_t)2},   {"s", (pdgid_t)3},
      {"s~", (pdgid_t)3},    {"c", (pdgid_t)4},     {"c~", (pdgid_t)4},   {"b", (pdgid_t)5},    {"b~", (pdgid_t)5},
      {"t", (pdgid_t)6},     {"t~", (pdgid_t)6},    {"e+", (pdgid_t)11},  {"e-", (pdgid_t)11},  {"ve", (pdgid_t)12},
      {"ve~", (pdgid_t)12},  {"mu+", (pdgid_t)13},  {"mu-", (pdgid_t)13}, {"vm", (pdgid_t)14},  {"vm~", (pdgid_t)14},
      {"tau+", (pdgid_t)15}, {"tau-", (pdgid_t)15}, {"vt", (pdgid_t)16},  {"vt~", (pdgid_t)16}, {"g", (pdgid_t)21},
      {"a", (pdgid_t)22},    {"z", (pdgid_t)23},    {"w+", (pdgid_t)24},  {"w-", (pdgid_t)24},  {"h", (pdgid_t)25},
  };

  MadGraphInterface::MadGraphInterface(const ParametersList& params)
      : SteeredObject(params),
        proc_(steer<std::string>("process")),
        model_(steer<std::string>("model")),
        tmp_dir_(steerAs<std::string, fs::path>("tmpDir")),
        card_path_(steerAs<std::string, fs::path>("cardPath")),
        log_filename_(steer<std::string>("logFile")),
        standalone_cpp_path_(steerAs<std::string, fs::path>("standaloneCppPath")),
        extra_particles_(steer<ParametersList>("extraParticles")) {
    if (proc_.empty() && standalone_cpp_path_.empty())
      throw CG_FATAL("MadGraphInterface") << "Neither a 'process' keyword nor a path to a MadGraph process interface "
                                             "already generated ('standaloneCppPath') was set to the parameters!\n"
                                          << params;
    std::ofstream log(log_filename_, std::ios::trunc);  // clearing the log
    parseExtraParticles();
  }

  void MadGraphInterface::parseExtraParticles() {
    for (const auto& extra_part : extra_particles_.keysOf<ParticleProperties>()) {
      if (mg5_parts_.count(extra_part) > 0)
        throw CG_FATAL("MadGraphInterface")
            << "Particle with name '" << extra_part << "' is already defined in internal LUT.";
      const auto& extra_part_prop = extra_particles_.get<ParticleProperties>(extra_part);
      // find the equivalent MadGraph particle to alias
      std::string found_mg_equiv;
      for (const auto& part : mg5_parts_)
        if (part.second == extra_part_prop.pdgid)
          found_mg_equiv = part.first;
      if (found_mg_equiv.empty())
        throw CG_FATAL("MadGraphInterface")
            << "No equivalent for particle with PDG id=" << extra_part_prop.pdgid << " in MadGraph LUT.";
      if (found_mg_equiv.at(found_mg_equiv.size() - 1) == '+' || found_mg_equiv.at(found_mg_equiv.size() - 1) == '-')
        found_mg_equiv.pop_back();
      if (extra_part_prop.charge != 0) {
        extra_part_definitions_ += "\ndefine " + extra_part + "+ = " + found_mg_equiv + "+";
        extra_part_definitions_ += "\ndefine " + extra_part + "- = " + found_mg_equiv + "-";
        mg5_parts_[extra_part + "+"] = mg5_parts_.at(found_mg_equiv + "+");
        mg5_parts_[extra_part + "-"] = mg5_parts_.at(found_mg_equiv + "-");
      } else {
        extra_part_definitions_ += "\ndefine " + extra_part + " = " + found_mg_equiv;
        mg5_parts_[extra_part] = mg5_parts_.at(found_mg_equiv);
      }
      //FIXME add extra particles properties (masses, ...)
      CG_DEBUG("MadGraphInterface") << "Defined '" << extra_part
                                    << "' as MadGraph alias for CepGen particle with properties: " << extra_part_prop
                                    << ".";
    }
    CG_LOG << extra_part_definitions_;
  }

  std::string MadGraphInterface::run() const {
    std::ofstream log(log_filename_, std::ios::app);  // appending at the end of the log

    fs::path cpp_path, cg_proc;
    if (!standalone_cpp_path_.empty()) {
      CG_INFO("MadGraphInterface:run") << "Running on a process already generated by mg5_aMC: " << standalone_cpp_path_;
      cpp_path = standalone_cpp_path_;
      cg_proc = tmp_dir_ / "cepgen_proc_interface.cpp";
    } else {
      CG_INFO("MadGraphInterface:run") << "Running the mg5_aMC process generation.";
      std::vector<std::string> cmds;
      if (!model_.empty()) {
        cmds.emplace_back("set auto_convert_model T");
        cmds.emplace_back("import model " + model_);
      }
      cmds.emplace_back(extra_part_definitions_);
      cmds.emplace_back("generate " + proc_);
      cmds.emplace_back("output standalone_cpp " + tmp_dir_.string());
      cpp_path = tmp_dir_;
      const auto num_removed_files = fs::remove_all(cpp_path);
      CG_DEBUG("MadGraphInterface:run") << "Removed " << utils::s("file", num_removed_files, true)
                                        << " from process directory " << cpp_path << ".";

      std::ofstream log(log_filename_, std::ios::app);  // appending at the end of the log
      log << "\n\n*** mg5_aMC process generation ***\n\n";
      log << utils::merge(mg5amc::runCommand(cmds, card_path_, true), "\n");

      CG_INFO("MadGraphInterface:run") << "Preparing the mg5_aMC process library.";
      log << "\n\n*** mg5_aMC process library compilation ***\n\n";
      cg_proc = prepareMadGraphProcess();
    }

#ifdef _WIN32
    fs::path lib_path = "CepGenMadGraphProcess.dll";
#else
    fs::path lib_path = "libCepGenMadGraphProcess.so";
#endif

    generateLibrary(cg_proc, cpp_path, lib_path);
    linkCards();
    return lib_path;
  }

  void MadGraphInterface::linkCards() const {
    for (const auto& f : fs::directory_iterator(tmp_dir_ / "Cards"))
      if (f.path().extension() == ".dat") {
        fs::path link_path = f.path().filename();
        if (!fs::exists(link_path))
          fs::create_symlink(f, link_path);
      }
    CG_DEBUG("MadGraphInterface:run") << "Created links in current directory for all cards in '" +
                                             std::string(tmp_dir_ / "Cards") + "'.";
  }

  std::string MadGraphInterface::prepareMadGraphProcess() const {
    //--- open template file
    std::ifstream tmpl_file(MADGRAPH_PROC_TMPL);
    std::string tmpl = std::string(std::istreambuf_iterator<char>(tmpl_file), std::istreambuf_iterator<char>());
    std::ofstream log(log_filename_, std::ios::app);  // appending at the end of the log
    log << "\n\n*** mg5_aMC process library compilation ***\n\n";

    const auto& parts = mg5amc::unpackProcessParticles(proc_);
    std::vector<pdgid_t> in_parts, out_parts;
    for (const auto& in_part : parts.first) {
      if (mg5_parts_.count(in_part) == 0) {
        const auto pprops = mg5amc::describeParticle(in_part, model_);
        mg5_parts_[in_part] = pprops.pdgid;
        PDG::get().define(pprops);
      }
      in_parts.emplace_back(mg5_parts_.at(in_part));
    }
    for (const auto& out_part : parts.second) {
      if (mg5_parts_.count(out_part) == 0) {
        const auto pprops = mg5amc::describeParticle(out_part, model_);
        mg5_parts_[out_part] = pprops.pdgid;
        PDG::get().define(pprops);
      }
      out_parts.emplace_back(mg5_parts_.at(out_part));
    }
    CG_INFO("MadGraphInterface.prepareMadGraphProcess") << "Unpacked process particles: "
                                                        << "incoming=" << in_parts << ", "
                                                        << "outgoing=" << out_parts << ".";

    const std::string process_description = proc_ + (!model_.empty() ? " (model: " + model_ + ")" : "");

    std::string src_filename = tmp_dir_ / "cepgen_proc_interface.cpp";
    std::ofstream src_file(src_filename);
    src_file << utils::replace_all(tmpl,
                                   {{"XXX_PART1_XXX", std::to_string(in_parts[0])},
                                    {"XXX_PART2_XXX", std::to_string(in_parts[1])},
                                    {"XXX_OUT_PART_XXX", utils::merge(out_parts, ", ")},
                                    {"XXX_PROC_NAME_XXX", mg5amc::normalise(proc_)},
                                    {"XXX_PROC_DESCRIPTION_XXX", process_description}});
    src_file.close();
    return src_filename;
  }

  void MadGraphInterface::generateLibrary(const fs::path& proc_path,
                                          const fs::path& in_path,
                                          const fs::path& out_lib) const {
    std::vector<std::string> src_files;
    src_files.emplace_back(proc_path.string());

    //--- find all processes registered
    std::vector<std::string> processes;
    try {
      for (const auto& p : fs::directory_iterator(in_path / "SubProcesses"))
        if (p.path().filename().string()[0] == 'P') {
          processes.emplace_back(p.path());
          for (const auto& f : fs::directory_iterator(p))
            if (f.path().extension() == ".cc")
              src_files.emplace_back(f.path());
        }
    } catch (const fs::filesystem_error& err) {
      throw CG_FATAL("MadGraphInterface:generateLibrary")
          << "Failed to retrieve all subprocesses in path " << in_path << "!\n"
          << err.what();
    }

    CG_DEBUG("MadGraphInterface:generateLibrary") << "Subprocess list: " << processes << ".";

    if (processes.size() != 1)
      throw CG_FATAL("MadGraphInterface:generateLibrary") << "Currently only single-process cases are supported!";

    //--- find all model source files
    for (const auto& f : fs::directory_iterator(in_path / "src"))
      if (f.path().extension() == ".cc")
        src_files.emplace_back(f.path());

#ifdef _WIN32
    throw CG_FATAL("MadGraphInterface:generateLibrary")
        << "Library generation not yet implemented for Window$ systems!";
#else
    utils::Caller::call({CC_CFLAGS,
                         "-fPIC",
                         "-shared",
                         "-Wno-unused-variable",
                         "-Wno-int-in-bool-context",
                         "-I" + (in_path / "src").string(),
                         "-I" + processes.at(0),
                         utils::merge(src_files, " "),
                         "-o " + out_lib.string()});
#endif
  }

  ParametersDescription MadGraphInterface::description() {
    auto desc = ParametersDescription();
    desc.add<std::string>("process", "").setDescription("MadGraph_aMC process definition");
    desc.add<std::string>("model", "sm-full").setDescription("MadGraph_aMC model name");
    desc.add<std::string>("cardPath", "/tmp/cepgen_mg5_input.dat")
        .setDescription("Temporary file where to store the input card for MadGraph_aMC");
    desc.add<std::string>("standaloneCppPath", "");
    desc.add<std::string>("tmpDir", "/tmp/cepgen_mg5_aMC")
        .setDescription("Temporary path where to store the MadGraph_aMC process definition files");
    desc.add<std::string>("logFile", "/tmp/cepgen_mg5_aMC.log")
        .setDescription("Temporary path where to store the log for this run");
    desc.add<ParametersDescription>("extraParticles", ParametersDescription())
        .setDescription("define internal MadGraph alias for a particle name");

    return desc;
  }
}  // namespace cepgen
