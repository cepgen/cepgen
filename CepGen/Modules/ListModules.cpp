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

#include <iomanip>

#include "CepGen/Core/Exception.h"
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
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  void dumpModules() {
    CG_LOG.log([](auto& info) {
      const std::string sep_mid(80, '-');
      info << "List of modules registered in the runtime database:";
      auto list_modules = [&info, &sep_mid](const auto& fact, const std::string& name) {
        info << "\n" << sep_mid << "\n" << utils::boldify(name);
        if (fact.empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : fact.modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << fact.describe(mod) << (fact.describeParameters(mod).empty() ? " (*)" : "");
      };
      auto list_int_modules = [&info, &sep_mid](const auto& fact,
                                                const std::string& name,
                                                std::function<std::string(int)> translator = nullptr) {
        info << "\n" << sep_mid << "\n" << utils::boldify(name);
        if (!translator)
          translator = [](int val) -> std::string { return std::to_string(val); };
        if (fact.empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : fact.modules())
          info << "\n> " << utils::colourise(translator(mod), utils::Colour::green, utils::Modifier::bold) << ": "
               << fact.describe(mod) << (fact.describeParameters(mod).empty() ? " (*)" : "");
      };

      list_modules(card::CardsHandlerFactory::get(), "Steering cards parsers");
      list_modules(IntegratorFactory::get(), "Integration algorithms");
      list_modules(proc::ProcessFactory::get(), "Physics processes");
      list_modules(formfac::FormFactorsFactory::get(), "Beam form factors modellings");
      list_int_modules(sigrat::SigmaRatiosFactory::get(), "Cross section ratios modellings");
      list_int_modules(strfun::StructureFunctionsFactory::get(), "Structure functions modellings", [](int mod) {
        std::ostringstream os;
        os << std::setw(3) << mod << "|" << (strfun::Type)mod;
        return os.str();
      });
      list_modules(EventModifierFactory::get(), "Event modification modules");
      list_modules(io::ExportModuleFactory::get(), "Export modules");
      list_modules(utils::FunctionalFactory::get(), "Functional evaluators");
      list_modules(utils::DrawerFactory::get(), "Drawer utilities");
      list_modules(AlphaEMFactory::get(), "alpha(EM) evolution algorithms");
      list_modules(AlphaSFactory::get(), "alpha(s) evolution algorithms");
    });
  }
}  // namespace cepgen
