/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGenAddOns/PythonWrapper/PythonConfigWriter.h"

namespace cepgen {
  namespace utils {
    PythonConfigWriter::PythonConfigWriter(const std::string& filename) : file_(filename) {
      file_ << "from sys import path\n"
            << "path.append('Cards')\n\n";
      file_ << "import Config.Core as cepgen\n\n";
    }

    PythonConfigWriter::~PythonConfigWriter() { file_.close(); }

    PythonConfigWriter& PythonConfigWriter::operator<<(const Parameters& params) {
      if (params.timeKeeper())
        (*this) << ParametersDescription("timer");
      if (params.hasProcess())
        (*this) << ParametersDescription(params.process().parameters()).setKey<std::string>("process");
      for (const auto& mod : params.eventModifiersSequence())
        (*this) << ParametersDescription(mod->parameters());
      for (const auto& mod : params.eventExportersSequence())
        (*this) << ParametersDescription(mod->parameters());
      return *this;
    }

    PythonConfigWriter& PythonConfigWriter::operator<<(const ParametersDescription& pdesc) {
      CG_DEBUG("PythonConfigWriter") << "Adding a parameters description object:\n" << pdesc;
      std::function<std::string(const ParametersDescription&, const std::string&, size_t)> write =
          [&](const ParametersDescription& pdesc, const std::string& key, size_t offset) -> std::string {
        // write a generic parameters description object
        std::stringstream os;
        const std::string off(offset * 4, ' ');
        os << off;
        if (!key.empty())
          os << key << " = ";

        const auto& params = pdesc.parameters();
        switch (pdesc.type()) {
          case ParametersDescription::Type::Module:
            os << "cepgen.Module(\"" << params.getNameString() << "\",";
            break;
          case ParametersDescription::Type::Parameters:
            os << "cepgen.Parameters(";
            break;
          case ParametersDescription::Type::ParametersVector:
            os << "list(";
            break;
          case ParametersDescription::Type::Value:
            break;
        }
        std::string sep = "";
        for (const auto& key : params.keys(false)) {
          os << sep << "\n";
          const auto& daugh = pdesc.get(key);
          switch (daugh.type()) {
            case ParametersDescription::Type::Module:
            case ParametersDescription::Type::Parameters:
              os << write(pdesc.get(key), key, offset + 1);
              break;
            case ParametersDescription::Type::ParametersVector: {
              std::string sep2;
              for (const auto& it : params.get<std::vector<ParametersList> >(key))
                os << sep2 << write(ParametersDescription(it), "", 0), sep2 = ", ";
            } break;
            case ParametersDescription::Type::Value:
              os << off << std::string(4, ' ') << key << " = ";
              if (params.has<bool>(key))
                os << (params.get<bool>(key) ? "True" : "False");
              else
                os << params.getString(key, true);
              break;
          }
          sep = ",";
        }
        switch (pdesc.type()) {
          case ParametersDescription::Type::Module:
          case ParametersDescription::Type::Parameters:
            os << "\n" << off;
            break;
          case ParametersDescription::Type::ParametersVector:
            os << ")" << off;
            break;
          case ParametersDescription::Type::Value:
            break;
        }
        os << ")";
        return os.str();
      };

      if (!pdesc.key().empty())
        file_ << pdesc.key() << " = ";
      file_ << write(pdesc, "", 0) << "\n";

      return *this;
    }
  }  // namespace utils
}  // namespace cepgen
