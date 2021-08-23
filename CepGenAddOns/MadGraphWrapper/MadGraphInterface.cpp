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
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphInterface.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

namespace cepgen {
  const std::unordered_map<std::string, pdgid_t> MadGraphInterface::mg5_parts_ = {
      {"d", (pdgid_t)1},     {"d~", (pdgid_t)1},    {"u", (pdgid_t)2},    {"u~", (pdgid_t)2},   {"s", (pdgid_t)3},
      {"s~", (pdgid_t)3},    {"c", (pdgid_t)4},     {"c~", (pdgid_t)4},   {"b", (pdgid_t)5},    {"b~", (pdgid_t)5},
      {"t", (pdgid_t)6},     {"t~", (pdgid_t)6},    {"e+", (pdgid_t)11},  {"e-", (pdgid_t)11},  {"ve", (pdgid_t)12},
      {"ve~", (pdgid_t)12},  {"mu+", (pdgid_t)13},  {"mu-", (pdgid_t)13}, {"vm", (pdgid_t)14},  {"vm~", (pdgid_t)14},
      {"tau+", (pdgid_t)15}, {"tau-", (pdgid_t)15}, {"vt", (pdgid_t)16},  {"vt~", (pdgid_t)16}, {"g", (pdgid_t)21},
      {"a", (pdgid_t)22},    {"z", (pdgid_t)23},    {"w+", (pdgid_t)24},  {"w-", (pdgid_t)24},  {"h", (pdgid_t)25},
  };

  MadGraphInterface::MadGraphInterface(const ParametersList& params)
      : proc_(params.get<std::string>("process")),
        model_(params.get<std::string>("model")),
        card_path_(params.get<std::string>("cardPath", "/tmp/cepgen_mg5_input.dat")),
        standalone_cpp_path_(params.get<std::string>("standaloneCppPath")),
        tmp_dir_(params.get<std::string>("tmpDir", "/tmp/cepgen_mg5_aMC")),
        log_filename_(params.get<std::string>("logFile", "/tmp/cepgen_mg5_aMC.log")) {
    if (proc_.empty())
      throw CG_FATAL("MadGraphInterface") << "'process' keyword not set to the parameters!\n" << params;
    std::ofstream log(log_filename_, std::ios::trunc);  // clearing the log
  }

  std::string MadGraphInterface::run() const {
    std::ofstream log(log_filename_, std::ios::app);  // appending at the end of the log

    std::string cpp_path;
    if (!standalone_cpp_path_.empty()) {
      CG_INFO("MadGraphInterface:run") << "Running on a process already generated by mg5_aMC:\n\t"
                                       << standalone_cpp_path_;
      cpp_path = standalone_cpp_path_;
    } else {
      CG_INFO("MadGraphInterface:run") << "Running the mg5_aMC process generation.";
      prepareCard();
      log << "\n\n*** mg5_aMC process generation ***\n\n";
      log << generateProcess(card_path_);
      cpp_path = tmp_dir_;
    }

#ifdef _WIN32
    std::string lib_path = "CepGenMadGraphProcess.dll";
#else
    std::string lib_path = "libCepGenMadGraphProcess.so";
#endif

    CG_INFO("MadGraphInterface:run") << "Preparing the mg5_aMC process library.";
    log << "\n\n*** mg5_aMC process library compilation ***\n\n";
    const auto cg_proc = prepareMadGraphProcess();
    log << generateLibrary(cg_proc, cpp_path, lib_path);

    CG_INFO("MadGraphInterface:run") << "Creating links for all cards in current directory.";
    linkCards();

    return lib_path;
  }

  void MadGraphInterface::prepareCard() const {
    std::ofstream card(card_path_);
    if (!model_.empty())
      card << "import model " << model_ << "\n";
    card << "generate " << proc_ << "\n";
    card << "output standalone_cpp " << tmp_dir_ << "\n";
    card << "exit\n";
    card.close();
  }

  void MadGraphInterface::linkCards() const {
    for (const auto& f : fs::directory_iterator(fs::path(tmp_dir_) / "Cards"))
      if (f.path().extension() == ".dat") {
        fs::path link_path = f.path().filename();
        if (!fs::exists(link_path))
          fs::create_symlink(f, link_path);
      }
  }

  std::string MadGraphInterface::prepareMadGraphProcess() const {
    //--- open template file
    std::ifstream tmpl_file(MADGRAPH_PROC_TMPL);
    std::string tmpl = std::string(std::istreambuf_iterator<char>(tmpl_file), std::istreambuf_iterator<char>());

    const auto& parts = unpackProcessParticles(proc_);

    const auto& in_parts = parts.first;
    utils::replace_all(tmpl, "XXX_PART1_XXX", std::to_string(in_parts[0]));
    utils::replace_all(tmpl, "XXX_PART2_XXX", std::to_string(in_parts[1]));

    const auto& out_parts = parts.second;
    std::ostringstream outparts_str;
    std::string sep;
    for (const auto& op : out_parts)
      outparts_str << sep << std::to_string(op), sep = ", ";
    utils::replace_all(tmpl, "XXX_OUT_PART_XXX", outparts_str.str());

    utils::replace_all(tmpl, "XXX_PROC_NAME_XXX", proc_);
    std::string descr = proc_;
    if (!model_.empty())
      descr += " (model: " + model_ + ")";
    utils::replace_all(tmpl, "XXX_PROC_DESCRIPTION_XXX", descr);

    std::string src_filename = fs::path(tmp_dir_) / "cepgen_proc_interface.cpp";
    std::ofstream src_file(src_filename);
    src_file << tmpl;
    src_file.close();
    return src_filename;
  }

  //--------------- static utilities ---------------

  MadGraphInterface::ProcessParticles MadGraphInterface::unpackProcessParticles(const std::string& proc) {
    ProcessParticles out;
    //--- dirty fix to specify incoming- and outgoing states
    //    as extracted from the mg5_aMC process string
    auto proc_name = proc;
    utils::trim(proc_name);
    const auto prim_proc = utils::split(proc_name, ',')[0];
    auto parts = utils::split(prim_proc, '>');
    if (parts.size() != 2)
      throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
          << "Unable to unpack particles from process name: \"" << proc_name << "\"";
    for (auto& p : parts)
      utils::trim(p);
    //--- incoming parton-like particles
    auto prim_parts = utils::split(parts[0], ' ');
    for (auto& p : prim_parts)
      p = utils::trim(p);
    CG_DEBUG("MadGraphInterface:unpackProcessParticles") << "Primary particles: " << prim_parts;
    if (prim_parts.size() != 2)
      throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
          << "Unable to unpack particles from process name: \"" << proc_name << "\"";
    for (const auto& p : prim_parts) {
      if (mg5_parts_.count(p) == 0)
        throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
            << "Particle with mg5_aMC name '" << p << "' was not recognised!";
      out.first.emplace_back(mg5_parts_.at(p));
    }
    //---- outgoing system
    auto dec_parts = utils::split(utils::trim(parts[1]), ' ');
    CG_DEBUG("MadGraphInterface:unpackProcessParticles") << "Outgoing system: " << dec_parts << ": " << parts;
    for (auto& p : dec_parts) {
      p = utils::trim(p);
      if (mg5_parts_.count(p) == 0)
        throw CG_FATAL("MadGraphInterface:unpackProcessParticles")
            << "Particle with mg5_aMC name '" << p << "' was not recognised!";
      out.second.emplace_back(mg5_parts_.at(p));
    }
    return out;
  }

  std::string MadGraphInterface::generateProcess(const std::string& in_path) {
    std::string cmd{MADGRAPH_BIN};
    cmd += " -f " + in_path;
    return runCommand(cmd);
  }

  std::string MadGraphInterface::generateLibrary(const std::string& proc_path,
                                                 const std::string& in_path,
                                                 const std::string& out_lib) {
    std::vector<std::string> src_files;
    src_files.emplace_back(proc_path);

    const fs::path tmp_path(in_path);

    //--- find all processes registered
    std::vector<std::string> processes;
    try {
      for (const auto& p : fs::directory_iterator(tmp_path / "SubProcesses"))
        if (p.path().filename().string()[0] == 'P') {
          processes.emplace_back(p.path());
          for (const auto& f : fs::directory_iterator(p))
            if (f.path().extension() == ".cc")
              src_files.emplace_back(f.path());
        }
    } catch (const fs::filesystem_error& err) {
      throw CG_FATAL("MadGraphInterface:generateLibrary")
          << "Failed to retrieve all subprocesses in path " << tmp_path << "!\n"
          << err.what();
    }

    CG_DEBUG("MadGraphInterface:generateLibrary") << "Subprocess list: " << processes << ".";

    if (processes.size() != 1)
      throw CG_FATAL("MadGraphInterface:generateLibrary") << "Currently only single-process cases are supported!";

    //--- find all model source files
    for (const auto& f : fs::directory_iterator(tmp_path / "src"))
      if (f.path().extension() == ".cc")
        src_files.emplace_back(f.path());

    std::string cmd{CC_CFLAGS};
    cmd += " -fPIC -shared";
    cmd += " -Wno-unused-variable -Wno-int-in-bool-context";
    cmd += " -I" + (tmp_path / "src").string();
    cmd += " -I" + processes.at(0);
    cmd += " " + utils::merge(src_files, " ");
    cmd += " -o " + out_lib;
    return runCommand(cmd);
  }

  std::string MadGraphInterface::runCommand(const std::string& cmd) {
    std::array<char, cmd_buffer_size_> buffer{};
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);

    CG_DEBUG("MadGraphInterface:runCommand") << "Running\n\t" << cmd;

    std::string result;
    while (fgets(buffer.data(), buffer.size(), pipe.get()))
      result += buffer.data();
    return result;
  }
}  // namespace cepgen
