/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#ifndef CepGen_Cards_LpairReader_h
#define CepGen_Cards_LpairReader_h

#include <memory>
#include <unordered_map>

#include "CepGen/Cards/Handler.h"

using std::string;

namespace cepgen {
  class ParametersList;
  namespace card {
    /// LPAIR-like steering cards parser and writer
    class LpairHandler final : public Handler {
    public:
      /// Read a LPAIR steering card
      explicit LpairHandler(const ParametersList&);
      static std::string description() { return "LPAIR-like cards parser"; }

      void pack(const Parameters*) override;
      Parameters* parse(const std::string&, Parameters*) override;
      /// Store a configuration into a LPAIR steering card
      void write(const std::string& file) const override;

    private:
      /// Single parameter handler
      /// \tparam T Parameter type
      template <typename T>
      struct Parameter {
        std::string key, description;
        T* value{nullptr};
      };
      /// Register a parameter to be steered to a configuration variable
      template <typename T>
      void registerParameter(const std::string& /*key*/, const std::string& /*description*/, T* /*def*/) {}
      template <typename T>
      void registerProcessParameter(const std::string& key,
                                    const std::string& description,
                                    const std::string& proc_key) {
        registerParameter<T>(key, description, &proc_params_->operator[]<T>(proc_key));
      }
      /// Register a kinematics block parameter to be steered
      template <typename T>
      void registerKinematicsParameter(const std::string& key,
                                       const std::string& description,
                                       const std::string& kin_key) {
        registerParameter<T>(key, description, &kin_params_->operator[]<T>(kin_key));
      }
      template <typename T>
      void registerGenerationParameter(const std::string& key,
                                       const std::string& description,
                                       const std::string& gen_key) {
        registerParameter<T>(key, description, &gen_params_->operator[]<T>(gen_key));
      }
      template <typename T>
      void registerIntegratorParameter(const std::string& key,
                                       const std::string& description,
                                       const std::string& int_key) {
        registerParameter<T>(key, description, &int_params_->operator[]<T>(int_key));
      }
      /// Set a parameter value
      template <typename T>
      void set(const std::string& /*key*/, const T& /*value*/) {}
      /// Retrieve a parameter value
      template <typename T>
      T get(const std::string& /*key*/) const {}

      void setParameter(const std::string& key, const std::string& value);
      std::string parameter(std::string key) const;
      std::string describe(std::string key) const;

      static constexpr int kInvalidInt = -999999;
      static constexpr double kInvalidDbl = 999.999;
      static constexpr const char* kInvalidStr = "(null)";

      std::unordered_map<std::string, Parameter<std::string> > p_strings_;
      std::unordered_map<std::string, Parameter<double> > p_doubles_;
      std::unordered_map<std::string, Parameter<int> > p_ints_;

      void init();
      std::shared_ptr<ParametersList> proc_params_, kin_params_, gen_params_, int_params_;
      int timer_{0}, iend_{1}, ext_log_{0};
      int str_fun_{11}, sr_type_{1}, lepton_id_{0};
      std::string proc_name_, evt_mod_name_, out_mod_name_;
      std::string out_file_name_, addons_list_;
      std::string kmr_grid_path_, mstw_grid_path_, pdg_input_path_;
    };

    //----- specialised registerers

    /// Register a string parameter
    template <>
    inline void LpairHandler::registerParameter<std::string>(const std::string& key,
                                                             const std::string& description,
                                                             std::string* def) {
      p_strings_[key] = Parameter<std::string>{key, description, def};
    }
    /// Register a double floating point parameter
    template <>
    inline void LpairHandler::registerParameter<double>(const std::string& key,
                                                        const std::string& description,
                                                        double* def) {
      p_doubles_[key] = Parameter<double>{key, description, def};
    }
    /// Register an integer parameter
    template <>
    inline void LpairHandler::registerParameter<int>(const std::string& key, const std::string& description, int* def) {
      p_ints_[key] = Parameter<int>{key, description, def};
    }

    //----- specialised setters

    template <>
    inline void LpairHandler::set<std::string>(const std::string& key, const std::string& value) {
      if (p_strings_.count(key))
        *p_strings_.at(key).value = value;
    }
    template <>
    inline void LpairHandler::set<double>(const std::string& key, const double& value) {
      if (p_doubles_.count(key))
        *p_doubles_.at(key).value = value;
    }
    template <>
    inline void LpairHandler::set<int>(const std::string& key, const int& value) {
      if (p_ints_.count(key))
        *p_ints_.at(key).value = value;
    }

    //----- specialised getters

    /// Retrieve a string parameter value
    template <>
    inline std::string LpairHandler::get(const std::string& key) const {
      if (p_strings_.count(key))
        return *p_strings_.at(key).value;
      return kInvalidStr;
    }
    /// Retrieve a floating point parameter value
    template <>
    inline double LpairHandler::get(const std::string& key) const {
      if (p_doubles_.count(key))
        return *p_doubles_.at(key).value;
      return -kInvalidDbl;
    }
    /// Retrieve an integer parameter value
    template <>
    inline int LpairHandler::get(const std::string& key) const {
      if (p_ints_.count(key))
        return *p_ints_.at(key).value;
      return -kInvalidInt;
    }
  }  // namespace card
}  // namespace cepgen

#endif
