/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

#include <fstream>

#include "CepGen/Generator.h"
#include "CepGen/Modules/DocumentationGeneratorFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/DocumentationGenerator.h"
#include "CepGen/Utils/Message.h"

using namespace std;

int main(int argc, char* argv[]) {
  string doc_generator, output_file;
  vector<string> categories, modules_names;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("documentation-generator,D", "type of documentation", &doc_generator, "text")
      .addOptionalArgument("output,o", "output file", &output_file, "")
      .addOptionalArgument("categories,C", "categories to document", &categories, vector<string>{})
      .addOptionalArgument("modules,m", "module names to document", &modules_names, vector<string>{})
      .parse();

  cepgen::initialise();
  auto gen = cepgen::DocumentationGeneratorFactory::get().build(cepgen::ParametersList()
                                                                    .setName(doc_generator)
                                                                    .feed(doc_generator)
                                                                    .set("categories", categories)
                                                                    .set("modules", modules_names));
  gen->initialise();
  const auto documentation = gen->describe();

  if (output_file.empty())
    CG_LOG << documentation;
  else {
    ofstream of(output_file);
    of << documentation;
    of.close();
    CG_INFO("main") << "Documentation written in '" << output_file << "'.";
  }

  return 0;
}
