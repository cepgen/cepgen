#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/ParametersDescription.h"
#include "CepGenAddOns/BoostWrapper/BoostTreeUtils.h"

namespace boost {
  namespace cepgen {
    pt::ptree pack(const ::cepgen::Parameters&) {
      pt::ptree out;
      return out;
    }

    pt::ptree pack(const ::cepgen::ParametersDescription& pdesc) { return pack(pdesc.parameters()); }

    pt::ptree pack(const ::cepgen::ParametersList& params) {
      pt::ptree out;
      for (const auto& key : params.keys()) {
        if (params.has<::cepgen::ParametersList>(key))
          out.add_child(key, pack(params.get<::cepgen::ParametersList>(key)));
        else if (params.has<int>(key))
          out.put(key, params.get<int>(key));
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
          throw ::cepgen::LoggedException(
              __FUNC__, "PythonConfigWriter", ::cepgen::Exception::Type::fatal, __FILE__, __LINE__)
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
        throw NullStream();
      for (const auto& it : tree) {
        try {
          if (it.first.empty())  // this might be a vector
            try {
              out.operator[]<std::vector<::cepgen::ParametersList>>(DAUGH_KEY).emplace_back(unpack(it.second));
            } catch (const boost::exception&) {
              try {
                out.operator[]<std::vector<double>>(DAUGH_KEY).emplace_back(it.second.get_value<double>());
              } catch (const boost::exception&) {
                try {
                  out.operator[]<std::vector<int>>(DAUGH_KEY).emplace_back(it.second.get_value<int>());
                } catch (const boost::exception&) {
                  out.operator[]<std::vector<std::string>>(DAUGH_KEY).emplace_back(it.second.get_value<std::string>());
                }
              }
            }
          else
            add(out, it.first, it.second);
        } catch (const Exception&) {
          if (it.second.get_value<std::string>().find('.') != std::string::npos)
            try {
              out.set<double>(it.first, it.second.get_value<double>());
            } catch (const boost::exception&) {
              out.set<std::string>(it.first, it.second.get_value<std::string>());
            }
          else
            try {
              out.set<int>(it.first, it.second.get_value<int>());
            } catch (const boost::exception&) {
              out.set<std::string>(it.first, it.second.get_value<std::string>());
            }
        }
      }
      return out;
    }
  }  // namespace cepgen
}  // namespace boost
