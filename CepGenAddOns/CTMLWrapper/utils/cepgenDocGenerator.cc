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

#include "CepGen/Generator.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/CTMLWrapper/DocumentationGenerator.h"

// list of factories documented
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

namespace cepgen {
  namespace utils {}
}  // namespace cepgen

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
  auto gen_params = cepgen::ParametersList()
                        .set<std::string>("output", output_file)
                        .set<bool>("useBS", use_bs)
                        .set<bool>("showGit", show_git)
                        .set<bool>("bare", bare);
  if (!show_title)
    gen_params.set<bool>("pageTitle", "");
  cepgen::utils::DocumentationGenerator gen{gen_params};

  gen.document("proc", "Processes", cepgen::ProcessFactory::get())
      //.document("cards", "Cards handler", cepgen::CardsHandlerFactory::get())
      .document("formfac", "Form factors", cepgen::FormFactorsFactory::get())
      .document("strfun", "Structure functions", cepgen::StructureFunctionsFactory::get())
      .document(
          "sigrat", "Longitudinal/transverse cross section ratio parameterisations", cepgen::SigmaRatiosFactory::get())
      .document("psmap", "Phase space mapper", cepgen::PhaseSpaceGeneratorFactory::get())
      .document("collflux", "Collinear parton flux modelling", cepgen::CollinearFluxFactory::get())
      .document("ktflux", "KT-factorised parton flux modelling", cepgen::KTFluxFactory::get())
      .document("alphaem", "Electromagnetic coupling evolution", cepgen::AlphaEMFactory::get())
      .document("alphas", "Strong coupling evolution", cepgen::AlphaSFactory::get())
      .document("integr", "Integrator algorithms", cepgen::IntegratorFactory::get())
      .document("func", "Functional parsers", cepgen::FunctionalFactory::get())
      .document("rndgen", "Random number generators", cepgen::RandomGeneratorFactory::get())
      .document("drawer", "Drawing tools", cepgen::DrawerFactory::get())
      .document("evtgen", "Event generation algorithms", cepgen::GeneratorWorkerFactory::get())
      .document("evtimp", "Event import algorithms", cepgen::EventImporterFactory::get())
      .document("evtmod", "Event modification algorithms", cepgen::EventModifierFactory::get())
      .document("evtout", "Event export modules", cepgen::EventExporterFactory::get());

  return 0;
}
