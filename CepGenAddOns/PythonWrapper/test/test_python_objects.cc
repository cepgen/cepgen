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

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::python::Environment env;

  {
    string str = "Héhéhé, test @ ünıc0d€";
    auto py_str = cepgen::python::set(str);
    auto str_ret = cepgen::python::get<string>(py_str.get());
    if (str_ret != str) {
      CG_LOG << "Object recasted from python is not identical to original object:\n"
             << "Original: " << str << ",\n"
             << "Recasted: " << str_ret << ".";
      return -1;
    }
    CG_LOG << "String test passed.";
  }
  {
    auto plist = cepgen::ParametersList().set<int>("foo", 42).set<double>("bar", M_PI).set<std::string>("baz", "héhé");
    auto py_dict = cepgen::python::set(plist);
    auto plist_ret = cepgen::python::get<cepgen::ParametersList>(py_dict.get());
    if (plist_ret != plist) {
      CG_LOG << "Object recasted from python is not identical to original object:\n"
             << "Original: " << plist << ",\n"
             << "Recasted: " << plist_ret << ".";
      return -1;
    }
    CG_LOG << "Parameters list/dictionary test passed.";
  }

  return 0;
}
