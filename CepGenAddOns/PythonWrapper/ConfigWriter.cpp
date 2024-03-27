/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/PythonWrapper/ConfigWriter.h"

using namespace std::string_literals;

namespace cepgen {
  namespace python {
    static std::string repr(const ParametersList& params, const std::string& key) {
      if (params.has<bool>(key))
        return params.get<bool>(key) ? "True" : "False";
      else if (params.has<int>(key))
        return "int(" + std::to_string(params.get<int>(key)) + ")";
      else if (params.has<unsigned long long>(key))
        return "int(" + std::to_string(params.get<unsigned long long>(key)) + ")";
      else if (params.has<std::string>(key))
        return "'" + utils::replaceAll(params.get<std::string>(key), "'", "\\'") + "'";
      else if (params.has<Limits>(key)) {
        const auto lim = params.get<Limits>(key);
        return "("s + std::to_string(lim.min()) + "," + (lim.hasMax() ? std::to_string(lim.max()) : "") + ")";
      } else if (params.has<std::vector<int> >(key))
        return "["s + utils::repr(params.get<std::vector<int> >(key)) + "]";
      else if (params.has<std::vector<double> >(key))
        return "["s + utils::repr(params.get<std::vector<double> >(key)) + "]";
      else if (params.has<std::vector<ParametersList> >(key)) {
        std::string out{"["}, sep;
        for (const auto& param : params.get<std::vector<ParametersList> >(key)) {
          out += sep + "cepgen.Parameters(";
          for (const auto& key : param.keys())
            out += key + " = " + repr(param, key);
          out += ")";
          sep = ", ";
        }
        return out + "]";
      }
      return params.getString(key, true);
    }

    ConfigWriter::ConfigWriter(const ParametersList& params)
        : SteeredObject(params), tab_len_(steer<int>("tabLength")) {
      if (steer<bool>("importPath"))
        os_ << "from sys import path\n"
            << "path.append('python')\n\n";
      os_ << "import Config.Core as cepgen\n\n";
    }

    ConfigWriter::~ConfigWriter() {
      if (const auto filename = steer<std::string>("filename"); !filename.empty()) {
        std::ofstream of(filename);
        of << os_.str();
        of.close();
      }
    }

    ConfigWriter& ConfigWriter::operator<<(const RunParameters& params) {
      if (params.timeKeeper())
        (*this) << ParametersDescription("timer");
      if (params.hasProcess())
        (*this) << ParametersDescription(params.process().parameters()).setKey<std::string>("process");
      for (const auto& mod : params.eventModifiersSequence())
        (*this) << ParametersDescription(mod->parameters()).setKey<std::string>("eventSequence");
      for (const auto& mod : params.eventExportersSequence())
        (*this) << ParametersDescription(mod->parameters()).setKey<std::string>("output");
      return *this;
    }

    ConfigWriter& ConfigWriter::operator<<(const ParametersDescription& pdesc) {
      CG_DEBUG("ConfigWriter") << "Adding a parameters description object:\n" << pdesc;
      const std::function<std::string(const ParametersDescription&, const std::string&, size_t)> write =
          [&](const ParametersDescription& pdesc, const std::string& key, size_t offset_num) -> std::string {
        // write a generic parameters description object
        std::stringstream os;
        os << offset(offset_num);
        if (!key.empty())
          os << key << " = ";

        std::string sep = "";
        const auto& params = pdesc.parameters();
        switch (pdesc.type()) {
          case ParametersDescription::Type::Module:
            os << "cepgen.Module("
               << (params.hasName<std::string>() ? "'" + params.getNameString() + "'"
                                                 : utils::toString(params.name<int>()));
            sep = ",";
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
        for (const auto& key : params.keys(false)) {
          os << sep << "\n";
          const auto& daugh = pdesc.get(key);
          switch (daugh.type()) {
            case ParametersDescription::Type::Module:
            case ParametersDescription::Type::Parameters:
              os << write(pdesc.get(key), key, offset_num + 1);
              break;
            case ParametersDescription::Type::ParametersVector: {
              os << offset(offset_num + 1) << key << " = [\n";
              for (const auto& it : params.get<std::vector<ParametersList> >(key))
                os << write(ParametersDescription(it), "", offset_num + 2) << ",\n";
              os << offset(offset_num + 1) << "]";
            } break;
            case ParametersDescription::Type::Value: {
              if (params.has<ParametersList>(key))
                os << write(ParametersDescription(params.get<ParametersList>(key)), key, offset_num + 1);
              else
                os << offset(offset_num + 1) << key << " = " << repr(params, key);
            } break;
          }
          sep = ",";
        }
        switch (pdesc.type()) {
          case ParametersDescription::Type::Module:
            if (!params.keys(false).empty())
              os << "\n" << offset(offset_num);
            break;
          case ParametersDescription::Type::Parameters:
            os << "\n" << offset(offset_num);
            break;
          case ParametersDescription::Type::ParametersVector:
            os << ")" << offset(offset_num);
            break;
          case ParametersDescription::Type::Value:
            break;
        }
        os << ")";
        return os.str();
      };
      const auto key = steer<bool>("camelCaseModuleNames") ? utils::toCamelCase(pdesc.key()) : pdesc.key();
      os_ << write(pdesc, key, 0) << "\n";
      return *this;
    }

    std::string ConfigWriter::operator()() const { return os_.str(); }

    ParametersDescription ConfigWriter::description() {
      auto desc = ParametersDescription();
      desc.add<bool>("importPath", false).setDescription("prepare the Python environment with path?");
      desc.add<bool>("camelCaseModuleNames", false).setDescription("convert the module names to camel case?");
      desc.add<int>("tabLength", 4).setDescription("number of spaces for one tabulation");
      desc.add<std::string>("filename", "").setDescription("Python output filename");
      return desc;
    }
  }  // namespace python
}  // namespace cepgen
