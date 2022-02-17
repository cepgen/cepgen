// clang-format off
#include "CepGenAddOns/PythonWrapper/PythonError.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/PythonWrapper/PythonConfigWriter.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::Generator gen;
  gen.parametersRef().setProcess(cepgen::proc::ProcessFactory::get().build("lpair"));
  cepgen::utils::PythonConfigWriter py("py_cfg.py");
  py << gen.parametersRef();
  return 0;
}
