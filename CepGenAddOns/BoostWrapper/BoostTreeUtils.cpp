/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2021  Laurent Forthomme
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
#include "CepGen/Parameters.h"
#include "CepGenAddOns/BoostWrapper/BoostTreeUtils.h"

namespace cepgen {
  namespace boost {
    pt::ptree pack(const Parameters&) {
      pt::ptree out;
      //FIXME implement this
      return out;
    }

    pt::ptree pack(const ParametersDescription& pdesc) { return pack(pdesc.parameters()); }

    pt::ptree pack(const ParametersList& params) {
      pt::ptree out;
      for (const auto& key : params.keys()) {
        if (params.has<ParametersList>(key))
          out.add_child(key, pack(params.get<ParametersList>(key)));
        else if (params.has<bool>(key))
          out.put(key, params.get<bool>(key));
        else if (params.has<int>(key))
          out.put(key, params.get<int>(key));
        else if (params.has<double>(key)) {
          std::ostringstream os;  // ensure floating point storage
          os << std::scientific << params.get<double>(key);
          out.put(key, os.str());
        } else if (params.has<std::string>(key))
          out.put(key, params.get<std::string>(key));
        else if (params.has<Limits>(key))
          out.add_child(key, pack(params.get<Limits>(key)));
        else if (params.has<std::vector<ParametersList>>(key))
          out.add_child(key, pack(params.get<std::vector<ParametersList>>(key)));
        else if (params.has<std::vector<int>>(key))
          out.add_child(key, pack(params.get<std::vector<int>>(key)));
        else if (params.has<std::vector<double>>(key))
          out.add_child(key, pack(params.get<std::vector<double>>(key)));
        else if (params.has<std::vector<std::string>>(key))
          out.add_child(key, pack(params.get<std::vector<std::string>>(key)));
        else
          throw CG_FATAL("BoostConfigWriter") << "Failed to recast the key \"" << key << "\" "
                                              << "with value \"" << params.getString(key) << "\"!";
      }
      return out;
    }

    template <>
    pt::ptree pack<ParametersList>(const std::vector<ParametersList>& vec) {
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

    pt::ptree pack(const Limits& lim) {
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

    ParametersList unpack(const pt::ptree& tree) {
      ParametersList out;
      if (tree.empty())
        return out;
      for (const auto& it : tree) {
        if (it.first.empty())  // this might be a vector
          try {
            out.operator[]<std::vector<ParametersList>>(DAUGH_KEY).emplace_back(unpack(it.second));
          } catch (const ::boost::exception&) {
            try {
              out.operator[]<std::vector<double>>(DAUGH_KEY).emplace_back(it.second.get_value<double>());
            } catch (const ::boost::exception&) {
              try {
                out.operator[]<std::vector<int>>(DAUGH_KEY).emplace_back(it.second.get_value<int>());
              } catch (const ::boost::exception&) {
                out.operator[]<std::vector<std::string>>(DAUGH_KEY).emplace_back(it.second.get_value<std::string>());
              }
            }
          }
        else if (it.second.get_value<std::string>().find('.') !=
                 std::string::npos)  // if contains a '.', might be a floating point variable
          try {
            out.set<double>(it.first, it.second.get_value<double>());
          } catch (const ::boost::exception&) {
            out.set<std::string>(it.first, it.second.get_value<std::string>());
          }
        else
          try {
            out.set<int>(it.first, it.second.get_value<int>());
          } catch (const ::boost::exception&) {
            out.set<std::string>(it.first, it.second.get_value<std::string>());
          }
      }
      CG_DEBUG("BoostTreeUtils:unpack") << "Unpacked parameters list:\n" << ParametersDescription(out) << ".";
      return out;
    }
  }  // namespace boost
}  // namespace cepgen
