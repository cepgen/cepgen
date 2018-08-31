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
#include "CepGen/Modules/AnalyticIntegratorFactory.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/CollinearFluxFactory.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/DerivatorFactory.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"

using namespace std;

#define LOOP_FACTORY(desc, obj)                                                                                 \
  if (all)                                                                                                      \
    os << "\n"                                                                                                  \
       << cepgen::utils::colourise(string(80, '=') + "\n" + string(desc) + " modules" + "\n" + string(80, '='), \
                                   cepgen::utils::Colour::green,                                                \
                                   cepgen::utils::Modifier::bold)                                               \
       << "\n";                                                                                                 \
  for (const auto& mod_name : obj::get().modules()) {                                                           \
    if (all)                                                                                                    \
      os << describe_one(desc, false, obj::get().describeParameters(mod_name));                                 \
    else                                                                                                        \
      for (const auto& mod : modules)                                                                           \
        if (mod == mod_name)                                                                                    \
          os << describe_one(desc, true, obj::get().describeParameters(mod_name));                              \
  }
#define LOOP_FACTORY_INT(desc, obj)                                                                             \
  if (all)                                                                                                      \
    os << "\n"                                                                                                  \
       << cepgen::utils::colourise(string(80, '=') + "\n" + string(desc) + " modules" + "\n" + string(80, '='), \
                                   cepgen::utils::Colour::green,                                                \
                                   cepgen::utils::Modifier::bold)                                               \
       << "\n";                                                                                                 \
  for (const auto& mod_name : obj::get().modules()) {                                                           \
    if (all)                                                                                                    \
      os << describe_one(desc, false, obj::get().describeParameters(mod_name));                                 \
    else                                                                                                        \
      for (const auto& mod : modules)                                                                           \
        if (cepgen::utils::isInt(mod) && stod(mod) == mod_name)                                                 \
          os << describe_one(desc, true, obj::get().describeParameters(mod_name));                              \
  }

/** Listing module for CepGen
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main(int argc, char* argv[]) {
  bool list_mods, safe_mode, dump_params, all;
  vector<string> addons, modules;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addOptionalArgument("list-modules,l", "list all runtime modules", &list_mods, false)
      .addOptionalArgument("modules,m", "list of runtime modules to be described", &modules)
      .addOptionalArgument("add-ons,e", "external runtime plugin", &addons)
      .addOptionalArgument("safe-mode,s", "safe mode", &safe_mode, false)
      .addOptionalArgument("dump-params,p", "dump the ParametersList object", &dump_params, false)
      .addOptionalArgument("all,a", "dump all modules descriptions", &all, false)
      .parse();

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
  if (all or !modules.empty()) {
    auto describe_one =
        [&dump_params](const string& type, bool dump_mod_name, const cepgen::ParametersDescription& desc) -> string {
      ostringstream os;
      os << "\n";
      if (dump_mod_name)
        os << cepgen::utils::colourise(
                  type, cepgen::utils::Colour::none, cepgen::utils::Modifier::underline | cepgen::utils::Modifier::bold)
           << " module:\n\n";
      os << desc.describe();
      if (dump_params)
        os << "\n\tParametersList object:\n\t\t" << desc.parameters();
      os << "\n";
      return os.str();
    };

    ostringstream os;
    LOOP_FACTORY("Cards steering", cepgen::CardsHandlerFactory)
    LOOP_FACTORY("Integrator", cepgen::IntegratorFactory)
    LOOP_FACTORY("Analytic integrator", cepgen::AnalyticIntegratorFactory)
    LOOP_FACTORY("Derivator", cepgen::DerivatorFactory)
    LOOP_FACTORY("Process", cepgen::ProcessFactory)
    LOOP_FACTORY("Parton flux modelling", cepgen::PartonFluxFactory)
    LOOP_FACTORY("Beam form factors modelling", cepgen::FormFactorsFactory)
    LOOP_FACTORY_INT("Collinear fluxes modelling", cepgen::CollinearFluxFactory)
    LOOP_FACTORY_INT("Structure functions modelling", cepgen::StructureFunctionsFactory)
    LOOP_FACTORY_INT("Cross sections ratio modelling", cepgen::SigmaRatiosFactory)
    LOOP_FACTORY("Event modification", cepgen::EventModifierFactory)
    LOOP_FACTORY("Export", cepgen::EventExporterFactory)
    LOOP_FACTORY("alpha(EM)", cepgen::AlphaEMFactory)
    LOOP_FACTORY("alpha(S)", cepgen::AlphaSFactory)
    LOOP_FACTORY("Functional evaluator", cepgen::FunctionalFactory)
    LOOP_FACTORY("Drawer", cepgen::DrawerFactory)
    CG_LOG << os.str();
    return 0;
  }

  return 0;
}
