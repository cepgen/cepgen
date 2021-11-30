#ifndef CepGen_Utils_ParametersDescription_h
#define CepGen_Utils_ParametersDescription_h

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  /// A description object for parameters collection
  class ParametersDescription : private ParametersList {
  public:
    /// Build the description of a parameters collection object
    /// \param[in] mod_name Module name (where applicable)
    explicit ParametersDescription(const std::string& mod_name = "");
    /// Build the (empty) description of a parameters collection object from its definition
    explicit ParametersDescription(const ParametersList& params);
    /// Does a description of this parameter (or parameters collection) exist?
    bool empty() const;
    /// Concatenate another description to this one
    ParametersDescription& operator+=(const ParametersDescription&);
    /// Set the description of this parameter (or parameters collection)
    ParametersDescription& setDescription(const std::string& descr);
    /// Description of this parameter (or parameters collection)
    const std::string& description() const { return mod_descr_; }
    /// Add the description to a new parameter
    template <typename T>
    ParametersDescription& add(const std::string& name, const T& def) {
      obj_descr_[name] = ParametersDescription();
      ParametersList::set<T>(name, def);
      return obj_descr_[name];
    }
    /// Human-readable description of all parameters and their default value
    std::string describe(size_t offset = 0) const;
    /// List of parameters associated to this description object
    const ParametersList& parameters() const { return *this; }

  private:
    std::string mod_descr_;
    std::map<std::string, ParametersDescription> obj_descr_;
  };
  /// Add the description to a new sub-description (aka ParametersList) object
  template <>
  ParametersDescription& ParametersDescription::add(const std::string&, const ParametersDescription&);
  /// Disable the addition of a ParametersList object to this description
  template <>
  ParametersDescription& ParametersDescription::add(const std::string&, const ParametersList&);
}  // namespace cepgen

#endif
