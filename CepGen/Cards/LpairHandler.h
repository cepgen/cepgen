#ifndef CepGen_Cards_LpairReader_h
#define CepGen_Cards_LpairReader_h

#include "CepGen/Cards/Handler.h"
#include <unordered_map>

using std::string;

namespace cepgen
{
  class ParametersList;
  namespace card
  {
    /// LPAIR-like steering cards parser and writer
    class LpairHandler : public Handler
    {
      public:
        /// Read a LPAIR steering card
        explicit LpairHandler( const char* file );

        /// Store a configuration into a LPAIR steering card
        void store( const char* file );
        static std::vector<std::string> split( const std::string&, char );

      private:
        /// Single parameter handler
        /// \tparam T Parameter type
        template<typename T> struct Parameter
        {
          Parameter( const char* key, const char* descr, T* value ) : key( key ), description( descr ), value( value ) {}
          std::string key, description;
          T* value;
        };
        /// Register a parameter to be steered to a configuration variable
        template<typename T> void registerParameter( const char* key, const char* description, T* def ) {}
        /// Set a parameter value
        template<typename T> void setValue( const char* key, const T& value ) {}
        /// Retrieve a parameter value
        template<typename T> T getValue( const char* key ) const {}

        void setParameter( const std::string& key, const std::string& value );
        std::string parameter( std::string key ) const;
        std::string description( std::string key ) const;

        static const int kInvalid;

        std::unordered_map<std::string, Parameter<std::string> > p_strings_;
        std::unordered_map<std::string, Parameter<double> > p_doubles_;
        std::unordered_map<std::string, Parameter<int> > p_ints_;
        std::unordered_map<std::string, Parameter<bool> > p_bools_;

        void init();
        std::shared_ptr<ParametersList> proc_params_;
        int str_fun_, sr_type_;
        double xi_min_, xi_max_;
        std::string proc_name_, evt_mod_name_, out_mod_name_;
        std::string out_file_name_;
        std::string integr_type_;
        std::string kmr_grid_path_, mstw_grid_path_;
        std::pair<unsigned short,unsigned short> hi_1_, hi_2_;
    };

    //----- specialised registerers

    /// Register a string parameter
    template<> inline void LpairHandler::registerParameter<std::string>( const char* key, const char* description, std::string* def ) { p_strings_.insert( std::make_pair( key, Parameter<std::string>( key, description, def ) ) ); }
    /// Register a double floating point parameter
    template<> inline void LpairHandler::registerParameter<double>( const char* key, const char* description, double* def ) { p_doubles_.insert( std::make_pair( key, Parameter<double>( key, description, def ) ) ); }
    /// Register an integer parameter
    template<> inline void LpairHandler::registerParameter<int>( const char* key, const char* description, int* def ) { p_ints_.insert( std::make_pair( key, Parameter<int>( key, description, def ) ) ); }
    /// Register a boolean parameter
    template<> inline void LpairHandler::registerParameter<bool>( const char* key, const char* description, bool* def ) { p_bools_.insert( std::make_pair( key, Parameter<bool>( key, description, def ) ) ); }

    //----- specialised setters

    template<> inline void LpairHandler::setValue<std::string>( const char* key, const std::string& value ) {
      auto it = p_strings_.find( key );
      if ( it != p_strings_.end() ) *it->second.value = value;
    }
    template<> inline void LpairHandler::setValue<double>( const char* key, const double& value ) {
      auto it = p_doubles_.find( key );
      if ( it != p_doubles_.end() ) *it->second.value = value;
    }
    template<> inline void LpairHandler::setValue<int>( const char* key, const int& value ) {
      auto it = p_ints_.find( key );
      if ( it != p_ints_.end() ) *it->second.value = value;
    }
    template<> inline void LpairHandler::setValue<bool>( const char* key, const bool& value ) {
      auto it = p_bools_.find( key );
      if ( it != p_bools_.end() ) *it->second.value = value;
    }

    //----- specialised getters

    /// Retrieve a string parameter value
    template<> inline std::string LpairHandler::getValue( const char* key ) const {
      const auto& it = p_strings_.find( key );
      if ( it != p_strings_.end() ) return *it->second.value;
      return "null";
    }
    /// Retrieve a floating point parameter value
    template<> inline double LpairHandler::getValue( const char* key ) const {
      const auto& it = p_doubles_.find( key );
      if ( it != p_doubles_.end() ) return *it->second.value;
      return -999.;
    }
    /// Retrieve an integer parameter value
    template<> inline int LpairHandler::getValue( const char* key ) const {
      const auto& it = p_ints_.find( key );
      if ( it != p_ints_.end() ) return *it->second.value;
      return 999;
    }
    /// Retrieve a boolean parameter value
    template<> inline bool LpairHandler::getValue( const char* key ) const {
      const auto& it = p_bools_.find( key );
      if ( it != p_bools_.end() ) return *it->second.value;
      return true;
    }
  }
}

#endif
