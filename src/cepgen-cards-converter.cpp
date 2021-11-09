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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Filesystem.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_config, output_config;
  bool debug;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addArgument("input,i", "input configuration", &input_config)
      .addArgument("output,o", "output output", &output_config)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();

  if (debug)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;

  cepgen::initialise();

  try {
    auto params = cepgen::card::Handler::parse(input_config);
    cepgen::card::Handler::write(params, output_config);
    CG_LOG << "Successfully converted the \"" << cepgen::utils::fileExtension(input_config) << "\" card into a \""
           << cepgen::utils::fileExtension(output_config) << "\" card.\n\t"
           << "\"" << output_config << "\" file created.";

  } catch (const cepgen::Exception& e) {
    throw CG_FATAL("main") << "Failed to convert a \"" << cepgen::utils::fileExtension(input_config)
                           << "\" card into a \"" << cepgen::utils::fileExtension(output_config) << "\" card!\n"
                           << e.message();
  }

  return 0;
}
