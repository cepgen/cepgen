/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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
#include "CepGenPython/ConfigWriter.h"

using namespace cepgen;
using namespace cepgen::python;
using namespace std::string_literals;

namespace cepgen::python {
  static std::string repr(const ParametersList& params, const std::string& key) {
    if (params.has<bool>(key))
      return params.get<bool>(key) ? "True" : "False";
    if (params.has<int>(key))
      return "int(" + std::to_string(params.get<int>(key)) + ")";
    if (params.has<unsigned long long>(key))
      return "int(" + std::to_string(params.get<unsigned long long>(key)) + ")";
    if (params.has<std::string>(key))
      return "'" + utils::replaceAll(params.get<std::string>(key), "'", "\\'") + "'";
    if (params.has<Limits>(key)) {
      const auto lim = params.get<Limits>(key);
      return "("s + std::to_string(lim.min()) + "," + (lim.hasMax() ? std::to_string(lim.max()) : "") + ")";
    }
    if (params.has<std::vector<Limits> >(key)) {
      std::string out{"["}, sep;
      for (const auto& lim : params.get<std::vector<Limits> >(key))
        out += sep + "("s + std::to_string(lim.min()) + "," + (lim.hasMax() ? std::to_string(lim.max()) : "") + ")",
            sep = ", ";
      return out + "]";
    }
    if (params.has<std::vector<int> >(key))
      return "["s + utils::repr(params.get<std::vector<int> >(key), ", ") + "]";
    if (params.has<std::vector<double> >(key))
      return "["s + utils::repr(params.get<std::vector<double> >(key), ", ") + "]";
    if (params.has<std::vector<std::vector<double> > >(key)) {
      std::string out{"["}, sep;
      for (const auto& vec : params.get<std::vector<std::vector<double> > >(key))
        out += sep + utils::repr(vec, ", "), sep = ", ";
      return out + "]";
    }
    if (params.has<std::vector<std::string> >(key))
      return "["s + utils::repr(params.get<std::vector<std::string> >(key), ", ") + "]";
    if (params.has<ParametersList>(key)) {
      const auto plist = params.get<ParametersList>(key);
      return (plist.hasName() ? "cepgen.Module(\'" + plist.name() + "\'" : "cepgen.Parameters(") + repr(plist, key) +
             ")";
    }
    if (params.has<std::vector<ParametersList> >(key)) {
      std::string out{"["}, sep;
      for (const auto& param : params.get<std::vector<ParametersList> >(key)) {
        out += sep + "cepgen.Parameters(";
        for (const auto& pkey : param.keys())
          out += pkey + " = " + repr(param, pkey);
        out += ")";
        sep = ", ";
      }
      return out + "]";
    }
    return params.getString(key, true);
  }
}  // namespace cepgen::python

ConfigWriter::ConfigWriter(const ParametersList& params) : SteeredObject(params), tab_len_(steer<int>("tabLength")) {
  if (steer<bool>("importPath"))
    os_ << "from sys import path\n"
        << "path.append('python')\n\n";
  os_ << "import Config.Core as cepgen\n\n";
}

ConfigWriter::~ConfigWriter() {
  if (const auto filename = steer<std::string>("filename"); !filename.empty()) {
    std::ofstream of(filename);
    of << os_.str();
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
      [&](const ParametersDescription& w_pdesc, const std::string& key, size_t offset_num) -> std::string {
    // write a generic parameters description object
    std::stringstream os;
    os << offset(offset_num);
    if (!key.empty())
      os << key << " = ";

    std::string sep = "";
    const auto& params = w_pdesc.parameters();
    switch (w_pdesc.type()) {
      case ParametersDescription::Type::Module:
        os << "cepgen.Module('" + params.getNameString() + "'";
        sep = ",";
        break;
      case ParametersDescription::Type::Value:
      case ParametersDescription::Type::Parameters:
        os << "cepgen.Parameters(";
        break;
      case ParametersDescription::Type::ParametersVector:
        os << "list(";
        break;
    }
    for (const auto& pkey : params.keys(false)) {
      os << sep << "\n";
      const auto& daugh = w_pdesc.get(pkey);
      switch (daugh.type()) {
        case ParametersDescription::Type::Module:
        case ParametersDescription::Type::Parameters:
          os << write(w_pdesc.get(pkey), pkey, offset_num + 1);
          break;
        case ParametersDescription::Type::ParametersVector: {
          os << offset(offset_num + 1) << pkey << " = [\n";
          for (const auto& it : params.get<std::vector<ParametersList> >(pkey))
            os << write(ParametersDescription(it), "", offset_num + 2) << ",\n";
          os << offset(offset_num + 1) << "]";
        } break;
        case ParametersDescription::Type::Value: {
          if (params.has<ParametersList>(pkey))
            os << write(ParametersDescription(params.get<ParametersList>(pkey)), pkey, offset_num + 1);
          else
            os << offset(offset_num + 1) << pkey << " = " << repr(params, pkey);
        } break;
      }
      sep = ",";
    }
    switch (w_pdesc.type()) {
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
  desc.add("importPath", false).setDescription("prepare the Python environment with path?");
  desc.add("camelCaseModuleNames", false).setDescription("convert the module names to camel case?");
  desc.add("tabLength", 4).setDescription("number of spaces for one tabulation");
  desc.add("filename", ""s).setDescription("Python output filename");
  return desc;
}
