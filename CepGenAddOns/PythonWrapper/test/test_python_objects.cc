#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  auto plist = cepgen::ParametersList().set<int>("foo", 42).set<double>("bar", M_PI).set<std::string>("baz", "héhé");
  CG_DEBUG("") << "Parameters list object to be \"dictionary-fied\": " << plist << ".";

  auto py_dict = cepgen::python::set(plist);

  return 0;
}
