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

  ParametersDescription::ParametersDescription(const ParametersDescription& oth)
      : ParametersList(oth), mod_descr_(oth.mod_descr_), obj_descr_(oth.obj_descr_) {}

  bool ParametersDescription::empty() const { return obj_descr_.empty() && mod_descr_.empty(); }

  ParametersDescription& ParametersDescription::operator=(const ParametersDescription& oth) {
    ParametersList::operator=(oth);
    mod_descr_ = oth.mod_descr_;
    obj_descr_ = oth.obj_descr_;
    return *this;
  }

  ParametersDescription& ParametersDescription::operator+=(const ParametersDescription& oth) {
    obj_descr_.insert(oth.obj_descr_.begin(), oth.obj_descr_.end());
    ParametersList::operator+=(oth);
    return *this;
  }

  std::string ParametersDescription::describe(size_t offset) const {
    static auto sep = [](size_t offset) -> std::string { return std::string(offset, '\t'); };
    const auto& mod_name = ParametersList::getString(ParametersList::MODULE_NAME);
    const auto& keys = ParametersList::keys(false);
    std::ostringstream os;
    if (mod_name.empty() && !keys.empty())
      os << utils::colourise("Parameters", utils::Colour::cyan, utils::Modifier::bold) << " collection ";
    else if (!mod_name.empty())
      os << utils::colourise("Module", utils::Colour::cyan, utils::Modifier::bold) << " " << utils::boldify(mod_name)
         << " ";
    if (!mod_descr_.empty())
      os << utils::colourise(mod_descr_, utils::Colour::none, utils::Modifier::italic);
    if (!keys.empty())
      os << "\n" << sep(offset + 1) << "List of parameters:";
    for (const auto& key : keys) {
      os << "\n" << sep(offset + 1) << "- " << utils::colourise(key, utils::Colour::none, utils::Modifier::underline);
      if (obj_descr_.count(key) > 0) {
        const auto& obj = obj_descr_.at(key);
        if (!ParametersList::has<ParametersList>(key))
          os << " (default value: " << ParametersList::getString(key) << ")";
        const auto& descr = obj.describe(offset + 1);
        if (!utils::trim(descr).empty())
          os << " " << descr;
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
    CG_DEBUG("ParametersDescription:add").log([this, &name, &desc](auto& log) {
      log << "Added a new parameters collection \"" << name << "\" as: " << desc.describe();
      const auto& mod_name = this->getString(ParametersList::MODULE_NAME);
      if (!mod_name.empty())
        log << "\nto the object with name: " << mod_name;
      log << ".";
    });
    return obj_descr_[name];
  }

  template <>
  ParametersDescription& ParametersDescription::add<ParametersList>(const std::string& name, const ParametersList&) {
    throw CG_FATAL("ParametersDescription:add")
        << "Using a ParametersList object for the description of a collection of parameters is not allowed.\n"
        << "Please use a ParametersDescription object instead for the description of the '" << name << "' collection.";
  }

  ParametersDescription& ParametersDescription::addParametersDescriptionVector(const std::string& name,
                                                                               const ParametersDescription& desc) {
    obj_descr_[name] = desc;
    ParametersList::set<std::vector<ParametersList> >(name, {});
    return obj_descr_[name];
  }

  ParametersList& ParametersDescription::parameters() { return *this; }

  const ParametersList& ParametersDescription::parameters() const { return *this; }

  void ParametersDescription::validate(const ParametersList&) const {
    throw CG_FATAL("ParametersDescription:validate") << "Not yet implemented!";
  }
}  // namespace cepgen
