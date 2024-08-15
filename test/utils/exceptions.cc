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

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  bool verbose;
  cepgen::ArgumentsParser(argc, argv).addOptionalArgument("verbose", "verbose mode", &verbose, false).parse();
  CG_TEST_DEBUG(verbose);

  const std::string test_string = "Haha, ceci est un test à géométrie variable! ☺";  // try with a bit of Unicode too
  for (int type = cepgen::Exception::Type::undefined; type < cepgen::Exception::Type::fatal; ++type) {
    ostringstream type_name;
    type_name << "Type " << static_cast<cepgen::Exception::Type>(type);
    auto throw_except = [&type, &test_string]() {
      throw cepgen::Exception("Test", "", static_cast<cepgen::Exception::Type>(type)) << test_string;
    };
    CG_TEST_EXCEPT(throw_except, type_name.str());
  }
  CG_TEST_SUMMARY;
}
