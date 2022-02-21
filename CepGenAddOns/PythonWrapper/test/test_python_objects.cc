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
#include "CepGenAddOns/PythonWrapper/PythonUtils.h"
// clang-format on
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

#define TEST_TYPE(type, object)                                                      \
  {                                                                                  \
    auto py_obj = cepgen::python::set(object);                                       \
    auto ret = cepgen::python::get<type>(py_obj.get());                              \
    if (ret != object) {                                                             \
      CG_LOG << "Object recasted from python is not identical to original object:\n" \
             << "Original: " << object << ",\n"                                      \
             << "Recasted: " << ret << ".";                                          \
      return -1;                                                                     \
    }                                                                                \
    CG_LOG << "type test passed.";                                                   \
  }

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::python::Environment env;
  TEST_TYPE(string, string("Héhéhé, test @ ünıc0d€"))
  TEST_TYPE(cepgen::ParametersList,
            cepgen::ParametersList().set<int>("foo", 42).set<double>("bar", M_PI).set<std::string>("baz", "héhé"))

  return 0;
}
