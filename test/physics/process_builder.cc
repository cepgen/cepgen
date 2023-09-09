/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

int main(int argc, char* argv[]) {
  string proc_name;
  bool list;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("proc-name,p", "name of the process", &proc_name, "lpair")
      .addOptionalArgument("list,l", "list all processes", &list, false)
      .parse();

  cepgen::initialise();

  if (list) {
    CG_LOG.log([](auto& log) {
      log << "List of modules registered in the runtime database:";
      for (const auto& mod : cepgen::ProcessFactory::get().modules())
        log << "\n> " << cepgen::utils::boldify(mod);
    });
    CG_TEST_SUMMARY;
  }

  if (proc_name.empty())
    CG_TEST_SUMMARY;

  CG_LOG << "Will build a process named \"" << proc_name << "\".";

  auto proc = cepgen::ProcessFactory::get().build(proc_name, cepgen::ParametersList());
  //--- at this point, the process has been found
  CG_LOG.log([&proc](auto& log) {
    log << "Successfully built the process \"" << proc->name() << "\"!\n"
        << " *) description: " << proc->description().description() << "\n"
        << " *) has event? " << proc->hasEvent() << "\n";
    if (proc->hasEvent()) {  //--- dump a typical event content
      log << "    event content (invalid kinematics, only check the parentage):\n";
      proc->event().dump();
    }
  });

  if (proc_name == "lpair") {
    CG_TEST_EQUAL(proc->hasEvent(), true, "LPAIR has event");
    if (proc->hasEvent()) {
      CG_TEST_EQUAL(proc->event().particles().size(), 9, "LPAIR particles content");
      const auto& cs = proc->event()(cepgen::Particle::Role::CentralSystem);
      CG_TEST_EQUAL(cs.size(), 2, "LPAIR outgoing state");
      CG_TEST_EQUAL(cs.at(0).pdgId(), 13, "LPAIR first outgoing particle");
      CG_TEST_EQUAL(cs.at(1).pdgId(), 13, "LPAIR second outgoing particle");
    }
  }

  CG_TEST_SUMMARY;
}
