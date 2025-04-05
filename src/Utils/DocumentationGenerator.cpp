/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/DocumentationGenerator.h"
#include "CepGen/Utils/String.h"

// list of factories documented
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/DocumentationGeneratorFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/GeneratorWorkerFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/PartonsPhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/PhaseSpaceGeneratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

using namespace cepgen;
using namespace cepgen::utils;

DocumentationGenerator::DocumentationGenerator(const ParametersList& params) : NamedModule(params) {
  const auto categories = steer<std::vector<std::string> >("categories"),
             mod_names = steer<std::vector<std::string> >("modules");
  const auto add_category =
      [this, &categories, &mod_names](
          const std::string& name, const std::string& title, const std::string& description, auto& factory) -> void {
    if (!categories.empty() && !contains(categories, name))
      return;
    category_t cat;
    cat.name = name;
    cat.title = title;
    cat.description = description;
    for (const auto& mod : factory.modules())
      if (mod_names.empty() || contains(mod_names, utils::toString(mod)) ||
          contains(mod_names, utils::toCamelCase(utils::toString(mod)))) {
        cat.modules[mod] = factory.describeParameters(mod).setKey(mod);
        for (const auto& [index, module_name] : factory.indices())
          if (module_name == mod)
            cat.modules_indices[mod] = index;
      }
    categories_.emplace_back(std::make_pair(name, cat));
  };
  add_category("proc", "Processes", "", cepgen::ProcessFactory::get());
  add_category("cards", "Cards handler", "", cepgen::CardsHandlerFactory::get());
  add_category("formfac", "Form factors", "", cepgen::FormFactorsFactory::get());
  add_category("strfun", "Structure functions", "", cepgen::StructureFunctionsFactory::get());
  add_category(
      "sigrat", "Longitudinal/transverse cross section ratio parameterisations", "", cepgen::SigmaRatiosFactory::get());
  add_category("partmap", "Partons generation algorithm", "", cepgen::PartonsPhaseSpaceGeneratorFactory::get());
  add_category("psmap", "Phase space mapper", "", cepgen::PhaseSpaceGeneratorFactory::get());
  add_category("collflux", "Collinear parton flux modelling", "", cepgen::CollinearFluxFactory::get());
  add_category("ktflux", "KT-factorised parton flux modelling", "", cepgen::KTFluxFactory::get());
  add_category("alphaem", "Electromagnetic coupling evolution", "", cepgen::AlphaEMFactory::get());
  add_category("alphas", "Strong coupling evolution", "", cepgen::AlphaSFactory::get());
  add_category("anaintegr", "Analytic integrator algorithms", "", cepgen::AnalyticIntegratorFactory::get());
  add_category("integr", "Integrator algorithms", "", cepgen::IntegratorFactory::get());
  add_category("func", "Functional parsers", "", cepgen::FunctionalFactory::get());
  add_category("rndgen", "Random number generators", "", cepgen::RandomGeneratorFactory::get());
  add_category("drawer", "Drawing tools", "", cepgen::DrawerFactory::get());
  add_category("evtgen", "Event generation algorithms", "", cepgen::GeneratorWorkerFactory::get());
  add_category("evtimp", "Event import algorithms", "", cepgen::EventImporterFactory::get());
  add_category("evtmod", "Event modification algorithms", "", cepgen::EventModifierFactory::get());
  add_category("evtout", "Event export modules", "", cepgen::EventExporterFactory::get());
  add_category("docs", "Documentation generator modules", "", cepgen::DocumentationGeneratorFactory::get());
}

ParametersDescription DocumentationGenerator::description() {
  auto desc = ParametersDescription();
  return desc;
}
