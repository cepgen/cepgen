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
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"

using namespace std;

/** Listing module for CepGen
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main(int argc, char* argv[]) {
  bool list_mods, debug, safe_mode, dump_params;
  vector<string> addons, modules;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addOptionalArgument("list-modules,l", "list all runtime modules", &list_mods, false)
      .addOptionalArgument("modules,m", "list of runtime modules to be described", &modules)
      .addOptionalArgument("add-ons,a", "external runtime plugin", &addons)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .addOptionalArgument("safe-mode,s", "safe mode", &safe_mode, false)
      .addOptionalArgument("dump-params,p", "dump the ParametersList object", &dump_params, false)
      .parse();

  //--- handle any debugging flag
  if (debug)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;

  //--- first start by defining the generator object
  try {
    for (const auto& lib : addons)
      cepgen::loadLibrary(lib);
    cepgen::Generator gen(safe_mode);
  } catch (const cepgen::Exception& e) {
    e.dump();
  }

  //--- if modules listing is requested
  if (list_mods) {
    cepgen::dumpModules();
    return 0;
  }
  if (!modules.empty()) {
    for (const auto& mod : modules) {
      auto describe = [&mod, &dump_params](const string& type, const cepgen::ParametersDescription& desc) {
        CG_LOG.log([&type, &desc, &mod, &dump_params](auto& log) {
          log << type << " module '" << mod << "'"
              << (desc.empty() ? " has no standard description" : ":\n" + desc.describe());
          if (dump_params)
            log << "\n\tParametersList object:\n\t\t" << desc.parameters();
        });
      };
      if (cepgen::card::CardsHandlerFactory::get().has(mod))
        describe("Cards steering", cepgen::card::CardsHandlerFactory::get().describeParameters(mod));
      if (cepgen::IntegratorFactory::get().has(mod))
        describe("Integrator", cepgen::IntegratorFactory::get().describeParameters(mod));
      if (cepgen::proc::ProcessFactory::get().has(mod))
        describe("Process", cepgen::proc::ProcessFactory::get().describeParameters(mod));
      if (cepgen::formfac::FormFactorsFactory::get().has(mod))
        describe("Beam form factors modelling", cepgen::formfac::FormFactorsFactory::get().describeParameters(mod));
      if (cepgen::utils::isNumber(mod)) {
        if (cepgen::strfun::StructureFunctionsFactory::get().has(stod(mod)))
          describe("Structure functions modelling",
                   cepgen::strfun::StructureFunctionsFactory::get().describeParameters(stod(mod)));
        if (cepgen::sigrat::SigmaRatiosFactory::get().has(stod(mod)))
          describe("Cross sections ratio modelling",
                   cepgen::sigrat::SigmaRatiosFactory::get().describeParameters(stod(mod)));
      }
      if (cepgen::EventModifierFactory::get().has(mod))
        describe("Event modification", cepgen::EventModifierFactory::get().describeParameters(mod));
      if (cepgen::io::ExportModuleFactory::get().has(mod))
        describe("Export", cepgen::io::ExportModuleFactory::get().describeParameters(mod));
      if (cepgen::utils::FunctionalFactory::get().has(mod))
        describe("Functional evaluator", cepgen::utils::FunctionalFactory::get().describeParameters(mod));
      if (cepgen::AlphaEMFactory::get().has(mod))
        describe("alpha(EM)", cepgen::AlphaEMFactory::get().describeParameters(mod));
      if (cepgen::AlphaSFactory::get().has(mod))
        describe("alpha(S)", cepgen::AlphaSFactory::get().describeParameters(mod));
    }
    return 0;
  }

  return 0;
}
