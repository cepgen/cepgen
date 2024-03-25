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

int main(int argc, char* argv[]) {
  std::string output_file;
  bool use_bs, show_title, show_git, bare;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("output,o", "output HTML file", &output_file, "index.html")
      .addOptionalArgument("bootstrap,b", "use Bootstrap CDN to prettify the output?", &use_bs, true)
      .addOptionalArgument("show-title,t", "show the page title?", &show_title, true)
      .addOptionalArgument("show-git,g", "show the git hash/branch?", &show_git, false)
      .addOptionalArgument("bare,e", "generate a bare version (without document tags) of the output?", &bare, false)
      .parse();

  cepgen::initialise();
  auto gen_params =
      cepgen::ParametersList().set<bool>("useBS", use_bs).set<bool>("showGit", show_git).set<bool>("bare", bare);
  if (!show_title)
    gen_params.set<bool>("pageTitle", "");
  auto gen = cepgen::DocumentationGeneratorFactory::get().build("ctml", gen_params);
  gen->initialise();
  const auto documentation = gen->describe();

  if (output_file.empty())
    CG_LOG << documentation;
  else {
    std::ofstream of(output_file);
    of << documentation;
    of.close();
    CG_INFO("main") << "Documentation written in '" << output_file << "'.";
  }

  return 0;
}
