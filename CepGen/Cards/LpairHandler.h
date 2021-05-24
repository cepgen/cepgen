#ifndef CepGen_Cards_LpairReader_h
#define CepGen_Cards_LpairReader_h

#include "CepGen/Cards/Handler.h"

#include <unordered_map>
#include <memory>

using std::string;

namespace cepgen {
  class ParametersList;
  namespace card {
    /// LPAIR-like steering cards parser and writer
    class LpairHandler : public Handler {
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
        T* value;
      };
      /// Register a parameter to be steered to a configuration variable
      template <typename T>
      void registerParameter(const std::string& key, const std::string& description, T* def) {}
      /// Register a kinematics block parameter to be steered
      template <typename T>
      void registerKinematicsParameter(const std::string& key,
                                       const std::string& description,
                                       const std::string& kin_key) {
        registerParameter<T>(key, description, &kin_params_->operator[]<T>(kin_key));
      }
      /// Set a parameter value
      template <typename T>
      void setValue(const std::string& key, const T& value) {}
      /// Retrieve a parameter value
      template <typename T>
      T getValue(const std::string& key) const {}

      void setParameter(const std::string& key, const std::string& value);
      std::string parameter(std::string key) const;
      std::string describe(std::string key) const;

      static const int kInvalid;

      std::unordered_map<std::string, Parameter<std::string> > p_strings_;
      std::unordered_map<std::string, Parameter<double> > p_doubles_;
      std::unordered_map<std::string, Parameter<int> > p_ints_;

      void init();
      std::shared_ptr<ParametersList> proc_params_, kin_params_, gen_params_;
      int timer_;
      int str_fun_, sr_type_, lepton_id_;
      std::string proc_name_, evt_mod_name_, out_mod_name_;
      std::string out_file_name_, addons_list_;
      std::string kmr_grid_path_, mstw_grid_path_, pdg_input_path_;
      int iend_;
    };

    //----- specialised registerers

    /// Register a string parameter
    template <>
    inline void LpairHandler::registerParameter<std::string>(const std::string& key,
                                                             const std::string& description,
                                                             std::string* def) {
      p_strings_.insert(std::make_pair(key, Parameter<std::string>{key, description, def}));
    }
    /// Register a double floating point parameter
    template <>
    inline void LpairHandler::registerParameter<double>(const std::string& key,
                                                        const std::string& description,
                                                        double* def) {
      p_doubles_.insert(std::make_pair(key, Parameter<double>{key, description, def}));
    }
    /// Register an integer parameter
    template <>
    inline void LpairHandler::registerParameter<int>(const std::string& key, const std::string& description, int* def) {
      p_ints_.insert(std::make_pair(key, Parameter<int>{key, description, def}));
    }

    //----- specialised setters

    template <>
    inline void LpairHandler::setValue<std::string>(const std::string& key, const std::string& value) {
      auto it = p_strings_.find(key);
      if (it != p_strings_.end())
        *it->second.value = value;
    }
    template <>
    inline void LpairHandler::setValue<double>(const std::string& key, const double& value) {
      auto it = p_doubles_.find(key);
      if (it != p_doubles_.end())
        *it->second.value = value;
    }
    template <>
    inline void LpairHandler::setValue<int>(const std::string& key, const int& value) {
      auto it = p_ints_.find(key);
      if (it != p_ints_.end())
        *it->second.value = value;
    }

    //----- specialised getters

    /// Retrieve a string parameter value
    template <>
    inline std::string LpairHandler::getValue(const std::string& key) const {
      const auto& it = p_strings_.find(key);
      if (it != p_strings_.end())
        return *it->second.value;
      return "null";
    }
    /// Retrieve a floating point parameter value
    template <>
    inline double LpairHandler::getValue(const std::string& key) const {
      const auto& it = p_doubles_.find(key);
      if (it != p_doubles_.end())
        return *it->second.value;
      return -999.;
    }
    /// Retrieve an integer parameter value
    template <>
    inline int LpairHandler::getValue(const std::string& key) const {
      const auto& it = p_ints_.find(key);
      if (it != p_ints_.end())
        return *it->second.value;
      return -999999;
    }
  }  // namespace card
}  // namespace cepgen

#endif
