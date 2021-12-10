#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/ParametersDescription.h"
#include "CepGen/Utils/PythonConfigWriter.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::initialise();
  cepgen::utils::PythonConfigWriter py("py_cfg.py");
  {
    auto strfun = cepgen::strfun::StructureFunctionsFactory::get().build(11);
    py << strfun->description();
  }
  return 0;
}
