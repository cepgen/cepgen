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

#include <iostream>

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"

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
    cout << "List of modules registered in the runtime database:";
    for (const auto& mod : cepgen::proc::ProcessFactory::get().modules())
      cout << "\n> " << cepgen::utils::boldify(mod);
    cout << endl;
    return 0;
  }

  if (!proc_name.empty()) {
    cout << "Will build a process named \"" << proc_name << "\"." << endl;

    auto proc = cepgen::proc::ProcessFactory::get().build(proc_name, cepgen::ParametersList());
    //--- at this point, the process has been found
    std::cout << "Successfully built the process \"" << proc->name() << "\"!\n"
              << " *) description: " << proc->parametersDescription().description() << "\n"
              << " *) has event? " << proc->hasEvent() << "\n";
    if (proc->hasEvent()) {  //--- dump a typical event content
      std::cout << "    event content (invalid kinematics, only check the parentage):\n";
      proc->addEventContent();
      proc->event().dump();
    }
  }

  return 0;
}
