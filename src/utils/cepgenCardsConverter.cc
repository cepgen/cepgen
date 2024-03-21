/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Filesystem.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_config, output_config;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addArgument("input,i", "input configuration", &input_config)
      .addArgument("output,o", "output output", &output_config)
      .parse();

  cepgen::initialise();

  try {
    auto in_card = cepgen::CardsHandlerFactory::get().parseFile(input_config);
    auto out_card = cepgen::CardsHandlerFactory::get().buildFromFilename(output_config);
    out_card->pack(in_card->runParameters());
    out_card->write(output_config);
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
