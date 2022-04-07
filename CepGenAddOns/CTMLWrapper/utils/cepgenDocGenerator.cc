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
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/CTMLWrapper/DocumentationGenerator.h"

namespace cepgen {
  namespace utils {}
}  // namespace cepgen

int main(int argc, char* argv[]) {
  std::string output_file;
  bool use_bs, bare;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("output,o", "output HTML file", &output_file, "index.html")
      .addOptionalArgument("bootstrap,b", "use Bootstrap CDN to prettify the output?", &use_bs, true)
      .addOptionalArgument("bare,e", "generate a bare version (without document tags) of the output?", &bare, false)
      .parse();

  cepgen::initialise();
  cepgen::utils::DocumentationGenerator gen{cepgen::ParametersList()
                                                .set<std::string>("output", output_file)
                                                .set<bool>("useBS", use_bs)
                                                .set<bool>("bare", bare)};

  gen.document("proc", "Processes", cepgen::proc::ProcessFactory::get())
      //.document("cards", "Cards handler", cepgen::card::CardsHandlerFactory::get())
      .document("formfac", "Form factors", cepgen::formfac::FormFactorsFactory::get())
      .document("strfun", "Structure functions", cepgen::strfun::StructureFunctionsFactory::get())
      .document("sigrat",
                "Longitudinal/transverse cross section ratio parameterisations",
                cepgen::sigrat::SigmaRatiosFactory::get())
      .document("alphaem", "Electromagnetic coupling evolution", cepgen::AlphaEMFactory::get())
      .document("alphas", "Strong coupling evolution", cepgen::AlphaSFactory::get())
      .document("integr", "Integrator algorithms", cepgen::IntegratorFactory::get())
      .document("func", "Functional parsers", cepgen::utils::FunctionalFactory::get())
      .document("drawer", "Drawing tools", cepgen::utils::DrawerFactory::get())
      .document("evtmod", "Event modification algorithms", cepgen::EventModifierFactory::get())
      .document("evtout", "Event export modules", cepgen::io::ExportModuleFactory::get());

  return 0;
}
