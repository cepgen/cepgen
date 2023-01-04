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

#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Version.h"

using namespace std;

int main(int argc, char* argv[]) {
  string commands;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("cmd,c", "list of arguments to be parsed", &commands, "")
      .parse();

  auto plist = cepgen::ParametersList().feed(commands);
  CG_LOG << "CepGen v" << cepgen::version::tag << " command line arguments parser parsed the following arguments:\n"
         << cepgen::ParametersDescription(plist);

  return 0;
}
