/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#include <pthread.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGenBoost/BoostTreeUtils.h"

namespace boost::cepgen {
  pt::ptree pack(const ::cepgen::RunParameters&) {
    pt::ptree out;
    return out;
  }

  pt::ptree pack(const ::cepgen::ParametersDescription& parameters_description) {
    return pack(parameters_description.parameters());
  }

  pt::ptree pack(const ::cepgen::ParametersList& params) {
    pt::ptree out;
    for (const auto& key : params.keys()) {
      if (params.has<::cepgen::ParametersList>(key))
        out.add_child(key, pack(params.get<::cepgen::ParametersList>(key)));
      else if (params.has<bool>(key))
        out.put(key, params.get<bool>(key));
      else if (params.has<int>(key))
        out.put(key, params.get<int>(key));
      else if (params.has<unsigned long long>(key))
        out.put(key, params.get<unsigned long long>(key));
      else if (params.has<double>(key)) {
        std::ostringstream os;  // ensure floating point storage
        os << std::scientific << params.get<double>(key);
        out.put(key, os.str());
      } else if (params.has<std::string>(key))
        out.put(key, params.get<std::string>(key));
      else if (params.has<::cepgen::Limits>(key))
        out.add_child(key, pack(params.get<::cepgen::Limits>(key)));
      else if (params.has<std::vector<::cepgen::ParametersList>>(key))
        out.add_child(key, pack(params.get<std::vector<::cepgen::ParametersList>>(key)));
      else if (params.has<std::vector<int>>(key))
        out.add_child(key, pack(params.get<std::vector<int>>(key)));
      else if (params.has<std::vector<double>>(key))
        out.add_child(key, pack(params.get<std::vector<double>>(key)));
      else if (params.has<std::vector<std::string>>(key))
        out.add_child(key, pack(params.get<std::vector<std::string>>(key)));
      else
        throw ::cepgen::Exception(__FUNC__, "boost::cepgen::pack", ::cepgen::Exception::Type::fatal, __FILE__, __LINE__)
            << "Failed to recast the key \"" << key << "\" "
            << "with value \"" << params.getString(key) << "\"!";
    }
    return out;
  }

  template <>
  pt::ptree pack<::cepgen::ParametersList>(const std::vector<::cepgen::ParametersList>& vec) {
    pt::ptree out;
    std::transform(vec.begin(), vec.end(), std::back_inserter(out), [](const auto& elem) {
      return std::make_pair("", pack(elem));
    });
    return out;
  }

  template <typename T>
  pt::ptree pack(const std::vector<T>& vec) {
    pt::ptree out;
    for (const auto& elem : vec) {
      pt::ptree elem_tree;
      elem_tree.put("", elem);
      out.push_back(std::make_pair("", elem_tree));
    }
    return out;
  }

  template <>
  pt::ptree pack<double>(const std::vector<double>& vec) {
    pt::ptree out;
    for (const auto& elem : vec) {
      pt::ptree elem_tree;
      std::ostringstream os;  // ensure floating point storage
      os << std::scientific << elem;
      elem_tree.put("", elem);
      out.push_back(std::make_pair("", elem_tree));
    }
    return out;
  }

  pt::ptree pack(const ::cepgen::Limits& lim) {
    pt::ptree out;
    if (lim.hasMin()) {
      pt::ptree min;
      std::ostringstream os;  // ensure floating point storage
      os << std::scientific << lim.min();
      min.put("", os.str());
      out.push_back(std::make_pair(MIN_KEY, min));
    }
    if (lim.hasMax()) {
      pt::ptree max;
      std::ostringstream os;  // ensure floating point storage
      os << std::scientific << lim.max();
      max.put("", os.str());
      out.push_back(std::make_pair(MAX_KEY, max));
    }
    return out;
  }

  ::cepgen::ParametersList unpack(const pt::ptree& tree) {
    ::cepgen::ParametersList out;
    if (tree.empty())
      return out;
    for (const auto& it : tree) {
      try {  // this might be a vector
        if (it.first.empty())
          try {  // try to parse as a parameters list
            out.operator[]<std::vector<::cepgen::ParametersList>>(DAUGH_KEY).emplace_back(unpack(it.second));
          } catch (const exception&) {
            try {  // try to parse as a float
              out.operator[]<std::vector<double>>(DAUGH_KEY).emplace_back(it.second.get_value<double>());
            } catch (const exception&) {
              try {  // try to parse as an integer
                out.operator[]<std::vector<int>>(DAUGH_KEY).emplace_back(it.second.get_value<int>());
              } catch (const exception&) {  // does not work, must be a string
                out.operator[]<std::vector<std::string>>(DAUGH_KEY).emplace_back(it.second.get_value<std::string>());
              }
            }
          }
        else
          add(out, it.first, it.second);
      } catch (const ::cepgen::Exception&) {
        if (it.second.get_value<std::string>().find('.') != std::string::npos)
          try {  // try to parse float if contains a '.'
            out.set<double>(it.first, it.second.get_value<double>());
          } catch (const exception&) {  // does not work as a float, must be a string
            out.set<std::string>(it.first, it.second.get_value<std::string>());
          }
        else
          try {  //  try to parse integer
            out.set<int>(it.first, it.second.get_value<int>());
          } catch (const exception&) {  // does not work as an integer, must be a string
            out.set<std::string>(it.first, it.second.get_value<std::string>());
          }
      }
    }
    return out;
  }

  void add(::cepgen::ParametersList& base, const std::string& name, const pt::ptree& tree) {
    auto plist = unpack(tree);
    //--- first check if we have a limits set
    if (plist.keys().size() <= 2 && (plist.has<double>(MIN_KEY) || plist.has<double>(MAX_KEY))) {
      // we might have a (min, max) limit
      ::cepgen::Limits lim;
      plist.fill<double>(MIN_KEY, lim.min());
      plist.fill<double>(MAX_KEY, lim.max());
      base.set<::cepgen::Limits>(name, lim);
    }
    //--- then check if daughter is a vector; if true, skip one hierarchy level
    else if (plist.has<std::vector<int>>(DAUGH_KEY))
      base.set<std::vector<int>>(name, plist.get<std::vector<int>>(DAUGH_KEY));
    else if (plist.has<std::vector<double>>(DAUGH_KEY)) {
      auto vec = plist.get<std::vector<double>>(DAUGH_KEY);
      base.set<std::vector<double>>(name, vec);
    } else if (plist.has<std::vector<std::string>>(DAUGH_KEY))
      base.set<std::vector<std::string>>(name, plist.get<std::vector<std::string>>(DAUGH_KEY));
    else
      base.set<::cepgen::ParametersList>(name, plist);
  }
}  // namespace boost::cepgen
