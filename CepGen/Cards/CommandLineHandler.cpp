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

#include <fstream>

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/TimeKeeper.h"

namespace cepgen {
  namespace card {
    /// Command line parser
    class CommandLineHandler final : public Handler {
    public:
      /// Cast command line arguments into a configuration word
      explicit CommandLineHandler(const ParametersList&);
      static std::string description() { return "Command line configuration parser"; }

      Parameters* parse(const std::string&, Parameters*) override;

    private:
      typedef std::vector<std::string> Args;

      static const double INVALID;

      Args argv_;
    };

    const double CommandLineHandler::INVALID = -999.999;

    CommandLineHandler::CommandLineHandler(const ParametersList& params)
        : Handler(params), argv_(params.get<std::vector<std::string> >("args")) {
      if (!filename_.empty())
        parse(filename_, rt_params_);
    }

    Parameters* CommandLineHandler::parse(const std::string& filename, Parameters* params) {
      if (!filename.empty()) {
        std::ifstream file(filename);
        std::string line;
        while (getline(file, line))
          argv_.emplace_back(line);
        file.close();
      }

      ParametersList pars;
      for (const auto& arg : argv_)
        pars.feed(arg);
      CG_INFO("CommandLineHandler") << "Arguments list: " << argv_ << " unpacked to:\n\t" << pars << ".";

      rt_params_ = params;

      //----- timer definition
      if (pars.get<bool>("timer", false))
        rt_params_->setTimeKeeper(new utils::TimeKeeper);

      //----- logging definition
      if (pars.has<int>("logging"))
        utils::Logger::get().level = pars.getAs<int, cepgen::utils::Logger::Level>("logging");
      else if (pars.has<ParametersList>("logging")) {
        const auto& log = pars.get<ParametersList>("logging");
        if (log.has<int>("level"))
          utils::Logger::get().level = log.getAs<int, cepgen::utils::Logger::Level>("level");
        utils::Logger::get().setExtended(log.get<bool>("extended", false));
      }

      //----- process definition
      auto proc = pars.get<ParametersList>("process");
      if (!proc.empty()) {
        if (rt_params_->hasProcess())
          proc = ParametersList(rt_params_->process().parameters()) + proc;
        rt_params_->setProcess(proc::ProcessFactory::get().build(proc));
      }

      //----- phase space definition
      auto kin = pars.get<ParametersList>("kinematics")
                     .set<ParametersList>("structureFunctions", pars.get<ParametersList>("strfun"))
                     .set<std::string>("formFactors", pars.get<std::string>("formfac"));
      rt_params_->kinematics = Kinematics(rt_params_->kinematics.parameters() + kin);

      //----- integration
      pars.fill<ParametersList>("integrator", *rt_params_->integrator);

      //----- events generation
      const auto& gen = pars.get<ParametersList>("generation");
      rt_params_->generation().setMaxGen(gen.get<int>("ngen", rt_params_->generation().maxGen()));
      if (gen.has<int>("nthreads"))
        rt_params_->generation().setNumThreads(gen.get<int>("nthreads"));
      if (gen.has<int>("nprn"))
        rt_params_->generation().setPrintEvery(gen.get<int>("nprn"));
      if (gen.has<int>("seed"))
        rt_params_->integrator->set<int>("seed", gen.get<int>("seed"));

      //----- event modification modules
      const auto& mod = pars.get<ParametersList>("eventmod");
      if (!mod.keys(true).empty()) {
        rt_params_->addModifier(EventModifierFactory::get().build(mod));
        rt_params_->eventModifiersSequence().rbegin()->get()->init();
      }

      //----- output modules definition
      const auto& out = pars.get<ParametersList>("output");
      if (!out.keys(true).empty())
        rt_params_->addOutputModule(io::ExportModuleFactory::get().build(out));
      return rt_params_;
    }
  }  // namespace card
}  // namespace cepgen

REGISTER_CARD_HANDLER(gCommandLineHandler, CommandLineHandler)
