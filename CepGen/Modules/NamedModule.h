#ifndef CepGen_Module_NamedModule_h
#define CepGen_Module_NamedModule_h

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  /// Base runtime module object
  template <typename T = std::string>
  class NamedModule {
  public:
    /// Build a module from its steering parameters
    explicit NamedModule(const ParametersList& params) : params_(params), name_(params.name<T>()) {}
    virtual ~NamedModule() = default;

    /// Module unique name
    const T& name() const { return name_; }
    /// Module description
    static inline std::string description() { return "No description"; }
    /// Collection of default parameters steering the module initialisation
    static inline ParametersList defaultParameters() { return ParametersList(); }
    /// Module user-defined parameters
    inline const ParametersList& parameters() const { return params_; }

  protected:
    /// Set of parameters to steer this output module
    const ParametersList params_;
    /// Module unique name
    const T name_;
  };
}  // namespace cepgen

#endif
