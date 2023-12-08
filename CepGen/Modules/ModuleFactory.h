/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2023  Laurent Forthomme
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
    static ModuleFactory& get();
    /// Disabled copy constructor
    ModuleFactory(const ModuleFactory&) = delete;
    /// Default destructor
    virtual ~ModuleFactory() = default;
    /// Disabled assignment operator
    void operator=(const ModuleFactory&) = delete;

    /// Describe the modules factory
    const std::string& description() const { return description_; }

    /// Register a named module in the database
    /// \tparam U Class to register (inherited from T base class)
    template <typename U>
    void registerModule(const I& name, const ParametersList& def_params = ParametersList()) {
      static_assert(std::is_base_of<T, U>::value,
                    "\n\n  *** Failed to register an object with improper inheritance into the factory. ***\n");
      if (has(name)) {
        std::ostringstream oss;
        oss << "\n\n  *** " << description_ << " detected a duplicate module registration for index/name \"" << name
            << "\"! ***\n";
        throw std::invalid_argument(oss.str());
      }
      map_.insert(std::make_pair(name, &buildModule<U>));
      auto desc = U::description();
      if (!def_params.empty())
        desc.parameters() += def_params;
      desc.parameters().setName(name);
      params_map_[name] = desc;
    }
    /// Build one instance of a named module
    /// \param[in] name Module name to retrieve
    /// \param[in] params List of parameters to be invoked by the constructor
    std::unique_ptr<T> build(const I& name, const ParametersList& params = ParametersList()) const;
    /// Build one instance of a named module
    /// \param[in] params List of parameters to be invoked by the constructor
    std::unique_ptr<T> build(const ParametersList& params = ParametersList()) const;

    /// Constructor type for a module
    typedef std::unique_ptr<T> (*Builder)(const ParametersList&);

    /// Describe one named module
    std::string describe(const I& name) const;
    /// Describe the parameters of one named module
    /// \params[in] params Parameters (incl. the name) to steer the description
    ParametersDescription describeParameters(const ParametersList&) const;
    /// Describe the parameters of one named module
    /// \params[in] name Name of the module to describe
    /// \params[in] params Additional parameters to steer the description
    ParametersDescription describeParameters(const I& name, const ParametersList& params = ParametersList()) const;

    /// List of modules registred in the database
    std::vector<I> modules() const;
    /// Is the database empty?
    bool empty() const { return map_.empty(); }
    /// Number of modules registered in the database
    size_t size() const { return map_.size(); }
    /// Check if a named module is registered
    bool has(const I& name) const { return map_.count(name) > 0; }

  private:
    /// Construct a module with its parameters set
    template <typename U>
    static std::unique_ptr<T> buildModule(const ParametersList& params) {
      return std::unique_ptr<T>(new U(params));
    }
    /// Factory name
    const std::string description_;
    /// Database of modules handled by this instance
    std::unordered_map<I, Builder> map_;
    /// Database of default parameters associated to modules
    std::unordered_map<I, ParametersDescription> params_map_;
    /// An empty parameters description
    const ParametersDescription empty_params_desc_;

  protected:
    /// Hidden default constructor for singleton operations
    explicit ModuleFactory(const std::string& descr = "Unnamed factory") : description_(descr) {}
  };
}  // namespace cepgen

#endif
