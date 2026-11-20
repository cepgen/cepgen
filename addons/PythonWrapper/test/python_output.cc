/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2026  Laurent Forthomme
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
#include "CepGenPython/Environment.h"
#include "CepGenPython/Error.h"
#include "CepGenPython/Utils.h"
// clang-format on
#include <CepGen/Core/Exception.h>
#include <CepGen/Core/ParametersDescription.h>
#include <CepGen/Core/RunParameters.h>
#include <CepGen/Generator.h>
#include <CepGen/Modules/FormFactorsFactory.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Process/Process.h>
#include <CepGen/Utils/ArgumentsParser.h>
#include <CepGen/Utils/Test.h>

#include <fstream>
#include <iterator>

#include "CepGenPython/ConfigWriter.h"

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
    gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(process));
    cepgen::python::ConfigWriter py(cepgen::ParametersList().set("filename", output_file));
    py << gen.runParameters();
  }

  try {
    auto env = cepgen::python::Environment(cepgen::ParametersList{});
    const auto path = cepgen::python::pythonPath(output_file);
    env.setProgramName(path);
    auto obj = cepgen::python::ObjectPtr::importModule(path);
    CG_TEST(obj != nullptr, "Module import");
    if (!obj)
      return -1;
    //CG_LOG << cepgen::python::ObjectPtr(PyObject_GenericGetDict(obj.get(), nullptr)).value<cepgen::ParametersList>();
    auto proc = obj.attribute("process");
    CG_TEST(proc != nullptr, "'process' attribute retrieval");
    const auto proc_params = proc.value<cepgen::ParametersList>();
    CG_TEST_EQUAL(proc_params.name(), process, "Process name conservation");
  } catch (const cepgen::python::Error& err) {
    err.dump();
  }
  {  // a single collection steered in place of a vector of collections (e.g. one form factors module for both beams)
    const auto ff_desc = cepgen::FormFactorsFactory::get().describeParameters("StandardDipole");
    cepgen::ParametersDescription desc("test");
    desc.addParametersDescriptionVector("formFactors", ff_desc, {ff_desc.parameters(), ff_desc.parameters()});
    desc = desc.steer(cepgen::ParametersList().set("formFactors", ff_desc.parameters()));
    const std::string vec_output_file = "py_vec_cfg.py";
    bool written = true;
    try {
      cepgen::python::ConfigWriter py(cepgen::ParametersList().set("filename", vec_output_file));
      py << desc;
    } catch (const cepgen::Exception& exc) {
      exc.dump();
      written = false;
    }
    CG_TEST(written, "writing a single collection steered in place of a vector");
    std::ifstream vec_file(vec_output_file);
    const std::string content((std::istreambuf_iterator<char>(vec_file)), std::istreambuf_iterator<char>());
    CG_TEST(content.find("formFactors = cepgen.Module('StandardDipole'") != std::string::npos,
            "single collection written in place of a vector");
  }
  CG_TEST_SUMMARY;
}
