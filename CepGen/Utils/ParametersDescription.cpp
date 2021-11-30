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

  bool ParametersDescription::empty() const { return obj_descr_.empty() && mod_descr_.empty(); }

  ParametersDescription& ParametersDescription::operator+=(const ParametersDescription& oth) {
    obj_descr_.insert(oth.obj_descr_.begin(), oth.obj_descr_.end());
    ParametersList::operator+=(oth);
    return *this;
  }

  std::string ParametersDescription::describe(size_t offset) const {
    std::ostringstream os;
    static auto sep = [](size_t offset) -> std::string { return std::string(offset, '\t'); };
    const auto& mod_name = ParametersList::name<std::string>();
    if (!mod_name.empty())
      os << sep(offset) << "Name: " << utils::boldify(mod_name) << "\n";
    if (!mod_descr_.empty())
      os << sep(offset + 1) << utils::colourise(mod_descr_, utils::Colour::reset, utils::Modifier::italic) << "\n";
    for (const auto& key : ParametersList::keys(false)) {
      os << sep(offset + 1) << utils::colourise(key, utils::Colour::reset, utils::Modifier::underline);
      if (obj_descr_.count(key) > 0) {
        const auto& obj = obj_descr_.at(key);
        if (!ParametersList::has<ParametersList>(key))
          os << " (default value: " << ParametersList::getString(key) << ")";
        os << "\n" << obj.describe(offset + 1);
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
