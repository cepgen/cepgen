/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

// clang-format off
#include "CepGenAddOns/PythonWrapper/PythonError.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/PythonConfigWriter.h"

using namespace std;

int main(int argc, char* argv[]) {
  string output_file, process;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("output,o", "Output python file", &output_file, "py_cfg.py")
      .addOptionalArgument("process,p", "Process name", &process, "")
      .parse();

  cepgen::Generator gen;
  if (process.empty())
    process = *cepgen::ProcessFactory::get().modules().begin();

  {
    gen.parametersRef().setProcess(cepgen::ProcessFactory::get().build(process));
    cepgen::utils::PythonConfigWriter py(output_file);
    py << gen.parametersRef();
  }

  try {
    auto env = cepgen::python::Environment();
    const auto path = cepgen::python::pythonPath(output_file);
    env.setProgramName(path);
    auto obj = cepgen::python::importModule(path);
    if (!obj) {
      CG_LOG << "Failed to import the module.";
      return -1;
    }
    //CG_LOG << cepgen::python::get<cepgen::ParametersList>(PyObject_GenericGetDict(obj.get(), nullptr));
    auto proc = cepgen::python::getAttribute(obj, "process");
    if (!proc) {
      CG_LOG << "Failed to retrieve a 'process' attribute.";
      return -1;
    }
    const auto proc_params = cepgen::python::get<cepgen::ParametersList>(proc);
    if (proc_params.name<std::string>() != process) {
      CG_LOG << "Process name was not conserved in Python file output.";
      return -1;
    }
  } catch (const cepgen::python::Error& err) {
    err.dump();
  }
  CG_LOG << "All tests passed.";

  return 0;
}
