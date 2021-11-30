#include <sstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/ParametersDescription.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  ParametersDescription::ParametersDescription(const std::string& mod_name) {
    if (!mod_name.empty())
      setName(mod_name);
  }

  ParametersDescription::ParametersDescription(const ParametersList& params) : ParametersList(params) {
    for (const auto& key : ParametersList::keys())
      obj_descr_[key] = ParametersDescription();
  }

  std::string ParametersDescription::describe(size_t offset) const {
    std::ostringstream os;
    static auto sep = [](size_t offset) -> std::string { return std::string(offset, '\t'); };
    const auto& mod_name = ParametersList::name<std::string>();
    if (!mod_name.empty())
      os << sep(offset) << "Name: " << utils::boldify(mod_name) << "\n";
    for (const auto& key : ParametersList::keys(false)) {
      os << sep(offset + 1) << utils::colourise(key, utils::Colour::reset, utils::Modifier::underline);
      if (obj_descr_.count(key) > 0) {
        const auto& obj = obj_descr_.at(key);
        if (!ParametersList::has<ParametersList>(key))
          os << " (default value: " << ParametersList::getString(key) << ")";
        os << "\n";
        if (!obj.description().empty())
          os << sep(offset + 2) << utils::colourise(obj.description(), utils::Colour::reset, utils::Modifier::italic)
             << "\n";
        os << obj.describe(offset + 2);
      }
    }
    return os.str();
  }

  ParametersDescription& ParametersDescription::setDescription(const std::string& descr) {
    mod_descr_ = descr;
    return *this;
  }

  template <>
  ParametersDescription& ParametersDescription::add<ParametersDescription>(const std::string& name,
                                                                           const ParametersDescription& desc) {
    obj_descr_[name] = desc;
    ParametersList::set<ParametersList>(name, desc.parameters());
    return obj_descr_[name];
  }

  template <>
  ParametersDescription& ParametersDescription::add<ParametersList>(const std::string& name, const ParametersList&) {
    throw CG_FATAL("ParametersDescription:add")
        << "Using a ParametersList object for the description of a collection of parameters is not allowed.\n"
        << "Please use a ParametersDescription object instead for the description of the '" << name << "' collection.";
  }
}  // namespace cepgen
