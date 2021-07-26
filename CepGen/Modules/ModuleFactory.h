#ifndef CepGen_Modules_ModuleFactory_h
#define CepGen_Modules_ModuleFactory_h

#include <memory>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "CepGen/Modules/NamedModule.h"

#define BUILDERNM(obj) obj##Builder
#define DEFINE_FACTORY_STR(name, objtype, descr)             \
  struct name : public ModuleFactory<objtype, std::string> { \
    explicit name() : ModuleFactory(descr) {}                \
    static name& get() {                                     \
      static name instance;                                  \
      return instance;                                       \
    }                                                        \
  }
#define DEFINE_FACTORY_INT(name, objtype, descr)     \
  struct name : public ModuleFactory<objtype, int> { \
    explicit name() : ModuleFactory(descr) {}        \
    static name& get() {                             \
      static name instance;                          \
      return instance;                               \
    }                                                \
  }

namespace cepgen {
  /// A generic factory to build modules
  /// \tparam T Base class to build
  /// \tparam I Indexing variable type
  template <typename T, typename I = std::string>
  class ModuleFactory {
  public:
    /// Retrieve a unique instance of this factory
    static ModuleFactory& get() {
      static ModuleFactory<T, I> instance;
      return instance;
    }
    /// Default destructor
    virtual ~ModuleFactory() = default;
    /// Register a named module in the database
    /// \tparam U Class to register (inherited from T base class)
    template <typename U>
    void registerModule(const I& name, const ParametersList& def_params = ParametersList()) {
      static_assert(std::is_base_of<T, U>::value,
                    "\n\n  *** Failed to register an object with improper inheritance into the factory. ***\n");
      if (map_.count(name) > 0 || params_map_.count(name) > 0) {
        std::ostringstream oss;
        oss << "\n\n  *** " << description_ << " detected a duplicate module registration for index/name \"" << name
            << "\"! ***\n";
        throw std::invalid_argument(oss.str());
      }
      map_[name] = &build<U>;
      descr_map_[name] = U::description();
      params_map_[name] = def_params;
    }
    /// Build one instance of a named module
    /// \param[in] name Module name to retrieve
    /// \param[in] params List of parameters to be invoked by the constructor
    std::unique_ptr<T> build(const I& name, ParametersList params = ParametersList()) const {
      if (name == I() || map_.count(name) == 0) {
        std::ostringstream oss;
        oss << "\n\n  *** " << description_ << " failed to build a module with index/name \"" << name << "\"! ***\n";
        throw std::invalid_argument(oss.str());
      }
      params.setName<I>(name);
      if (params_map_.count(name) > 0)
        params += params_map_.at(name);
      return map_.at(name)(params);
    }
    /// Build one instance of a named module
    /// \param[in] params List of parameters to be invoked by the constructor
    std::unique_ptr<T> build(ParametersList params = ParametersList()) const {
      if (params.has<I>(ParametersList::MODULE_NAME)) {
        const I& idx = params.get<I>(ParametersList::MODULE_NAME);
        if (map_.count(idx) == 0) {
          std::ostringstream oss;
          oss << "\n\n  *** " << description_ << " failed to build a module with index/name \"" << idx << "\"! ***\n";
          throw std::invalid_argument(oss.str());
        }
        if (params_map_.count(idx) > 0)
          params += params_map_.at(idx);
        return map_.at(idx)(params);
      } else
        throw std::invalid_argument("\n\n  *** " + description_ +
                                    " failed to retrieve an indexing key from parameters to build the module! ***\n");
    }
    /// Describe the modules factory
    const std::string& description() { return description_; }
    /// Describe one named module
    const std::string& describe(const I& name) const { return descr_map_.at(name); }
    /// List of modules registred in the database
    std::vector<I> modules() const {
      std::vector<I> out;
      for (const auto& p : map_)
        out.emplace_back(p.first);
      return out;
    }

  private:
    /// Construct a module with its parameters set
    template <typename U>
    static std::unique_ptr<T> build(const ParametersList& params) {
      return std::unique_ptr<T>(new U(params));
    }
    /// Factory name
    const std::string description_;
    /// Constructor type for a module
    typedef std::unique_ptr<T> (*ModCreate)(const ParametersList&);
    /// Database of modules handled by this instance
    std::unordered_map<I, ModCreate> map_;
    /// Database of default-constructed objects
    std::unordered_map<I, std::string> descr_map_;
    /// Database of default parameters associated to modules
    std::unordered_map<I, ParametersList> params_map_;

  protected:
    explicit ModuleFactory(const std::string& descr = "Unnamed factory") : description_(descr) {}

  public:
    ModuleFactory(const ModuleFactory&) = delete;
    void operator=(const ModuleFactory&) = delete;
  };
}  // namespace cepgen

#endif
