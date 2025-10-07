/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::initialise();

  bool include_mg5_proc;
  vector<string> processes;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("processes,p", "name of the processes", &processes, cepgen::ProcessFactory::get().modules())
      .addOptionalArgument("include-madgraph"s, "include MG5_aMC process?", &include_mg5_proc, true)
      .parse();

  CG_LOG << "Will test process(es): " << processes << ".";
  for (const auto& proc_name : processes) {
    auto proc_name_fix = proc_name;
    if (const auto mg5_proc = cepgen::utils::split(proc_name, ':'); mg5_proc.at(0) == "mg5_aMC") {
      if (!include_mg5_proc)
        continue;
      if (mg5_proc.size() > 1)
        proc_name_fix += "<process:'a e- > mu- mu+ e-'<removeLibrary:true";
      else
        proc_name_fix += "<process:'a a > mu- mu+'<removeLibrary:true";
    }
    auto proc = cepgen::ProcessFactory::get().build(proc_name_fix);
    proc->initialise();
    CG_DEBUG("main").log([&proc](auto& log) {
      log << "Successfully built the process \"" << proc->name() << "\"!\n"
          << " *) description: " << proc->description().description() << "\n"
          << " *) has event? " << proc->hasEvent() << "\n";
      if (proc->hasEvent()) {  // dump a typical event content
        log << "    event content (invalid kinematics, only check the parentage):\n";
        proc->event().dump();
      }
    });
    CG_TEST_EQUAL(proc->hasEvent(), true, "process has event");
    if (!proc->hasEvent())
      continue;
    if (proc_name == "lpair"s || proc_name == "pptoff"s || proc_name == "mg5_aMC"s) {
      CG_TEST_EQUAL(proc->event().particles().size(), 9, proc_name + " particles content");
      const auto& cs = proc->event()(cepgen::Particle::Role::CentralSystem);
      CG_TEST_EQUAL(cs.size(), 2, proc_name + " outgoing state");
      CG_TEST_EQUAL(cs.at(0).integerPdgId(), +13, proc_name + " first outgoing particle");
      CG_TEST_EQUAL(cs.at(1).integerPdgId(), -13, proc_name + " second outgoing particle");
    }
    if (proc_name == "pptoww"s) {
      CG_TEST_EQUAL(proc->event().particles().size(), 9, proc_name + " particles content");
      const auto& cs = proc->event()(cepgen::Particle::Role::CentralSystem);
      CG_TEST_EQUAL(cs.size(), 2, proc_name + " outgoing state");
      CG_TEST_EQUAL(cs.at(0).integerPdgId(), +24, proc_name + " first outgoing particle");
      CG_TEST_EQUAL(cs.at(1).integerPdgId(), -24, proc_name + " second outgoing particle");
    }
  }

  CG_TEST_SUMMARY;
}
