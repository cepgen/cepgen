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
#include <sstream>

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  namespace io {
    ExportModule::ExportModule(const ParametersList& params) : NamedModule(params) {}

    ExportModule::~ExportModule() {
      CG_DEBUG("ExportModule") << "Destructor called for output module \"" << name_ << "\".";
    }

    std::string ExportModule::banner(const Parameters& params, const std::string& prep) {
      const size_t len = 45 + version::tag.size();
      std::ostringstream os;
      os << prep << "******* Sample generated with CepGen " << version::tag << " *******\n"
         << prep << " Process: " << params.processName() << " (" << params.kinematics().incomingBeams().mode() << ")\n";
      if (params.kinematics().incomingBeams().mode() != mode::Kinematics::ElasticElastic)
        os << prep << " Structure functions: "
           << params.kinematics().incomingBeams().structureFunctions()->parametersDescription().description() << "\n";
      if (!params.eventModifiersSequence().empty()) {
        os << prep << " " << utils::s("Event modifier", params.eventModifiersSequence().size()) << ": ";
        std::string sep;
        for (const auto& mod : params.eventModifiersSequence())
          os << sep << mod->name(), sep = ", ";
        os << "\n";
      }
      const auto& cuts = params.kinematics().cuts();
      os << prep << std::left << std::setw(len) << std::setfill('*') << "*** Incoming state "
         << "\n";
      for (const auto& cut : cuts.initial.list())
        os << prep << " " << cut.description << ": " << cut.limits << "\n";
      os << prep << std::setw(len) << std::setfill('*') << "*** Central system "
         << "\n";
      for (const auto& cut : cuts.central.list())
        os << prep << " " << cut.description << ": " << cut.limits << "\n";
      if (params.kinematics().incomingBeams().mode() != mode::Kinematics::ElasticElastic) {
        os << prep << std::setw(len) << std::setfill('*') << "*** Remnants states "
           << "\n";
        for (const auto& cut : cuts.remnants.list())
          os << prep << " " << cut.description << ": " << cut.limits << "\n";
      }
      os << prep << std::string(45 + version::tag.size(), '*');
      return os.str();
    }

    ParametersDescription ExportModule::parametersDescription() {
      auto desc = ParametersDescription();
      return desc;
    }
  }  // namespace io
}  // namespace cepgen
