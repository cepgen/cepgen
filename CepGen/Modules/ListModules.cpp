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

#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Modules/IntegratorFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/AlphaS.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  void dumpModules() {
    CG_LOG.log([](auto& info) {
      const std::string sep_mid(80, '-');
      info << "List of modules registered in the runtime database:\n";
      {
        info << sep_mid << "\n" << utils::boldify("Steering cards parsers");
        if (card::CardsHandlerFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : card::CardsHandlerFactory::get().modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold)
               << " extension: " << card::CardsHandlerFactory::get().describe(mod);
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Integration algorithms");
        if (IntegratorFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : IntegratorFactory::get().modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << IntegratorFactory::get().describe(mod);
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Physics processes");
        if (proc::ProcessFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : proc::ProcessFactory::get().modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << proc::ProcessFactory::get().describe(mod);
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Beam form factors modellings");
        if (formfac::FormFactorsFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : formfac::FormFactorsFactory::get().modules()) {
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << formfac::FormFactorsFactory::get().describe(mod);
        }
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Structure functions modellings");
        if (strfun::StructureFunctionsFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : strfun::StructureFunctionsFactory::get().modules()) {
          std::ostringstream os;
          os << (strfun::Type)mod;
          info << "\n> " << utils::colourise(std::to_string(mod), utils::Colour::green, utils::Modifier::bold) << "|"
               << utils::colourise(os.str(), utils::Colour::green, utils::Modifier::bold) << ": "
               << strfun::StructureFunctionsFactory::get().describe(mod);
        }
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Cross section ratios modellings");
        if (sigrat::SigmaRatiosFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : sigrat::SigmaRatiosFactory::get().modules())
          info << "\n> " << utils::colourise(std::to_string(mod), utils::Colour::green, utils::Modifier::bold) << ": "
               << sigrat::SigmaRatiosFactory::get().describe(mod);
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Event modification modules");
        if (EventModifierFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : EventModifierFactory::get().modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << EventModifierFactory::get().describe(mod);
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Export modules");
        if (io::ExportModuleFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : io::ExportModuleFactory::get().modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << io::ExportModuleFactory::get().describe(mod);
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("Functional evaluators");
        if (utils::FunctionalFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : utils::FunctionalFactory::get().modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << utils::FunctionalFactory::get().describe(mod);
      }
      {
        info << "\n" << sep_mid << "\n" << utils::boldify("alpha(s) evolution algorithms");
        if (AlphaSFactory::get().modules().empty())
          info << "\n>>> " << utils::colourise("none found", utils::Colour::red) << " <<<";
        for (const auto& mod : AlphaSFactory::get().modules())
          info << "\n> " << utils::colourise(mod, utils::Colour::green, utils::Modifier::bold) << ": "
               << AlphaSFactory::get().describe(mod);
      }
    });
  }
}  // namespace cepgen
