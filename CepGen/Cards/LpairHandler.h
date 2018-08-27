#ifndef CepGen_Cards_LpairReader_h
#define CepGen_Cards_LpairReader_h

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/ParametersList.h"
#include <unordered_map>

using std::string;

namespace CepGen
{
  namespace Cards
  {
    /// LPAIR-like steering cards parser and writer
    class LpairHandler : public Handler
    {
      public:
        /// Read a LPAIR steering card
        explicit LpairHandler( const char* file );

        /// Store a configuration into a LPAIR steering card
        void store( const char* file );

      private:
        template<class T> struct Parameter {
          Parameter( const char* key, const char* descr, T* value ) : key( key ), description( descr ), value( value ) {}
          std::string key, description;
          T* value;
        };
        /// Register a parameter to be steered to a configuration variable
        template<class T> void registerParameter( const char* key, const char* description, T* def ) {}
        /// Set a parameter value
        template<class T> void setValue( const char* key, const T& value ) {}
        /// Retrieve a parameter value
        template<class T> T getValue( const char* key ) const {}

        void setParameter( const std::string& key, const std::string& value );
        std::string getParameter( std::string key ) const;
        std::string getDescription( std::string key ) const;

        static const int kInvalid;

        std::unordered_map<std::string, Parameter<std::string> > p_strings_;
        std::unordered_map<std::string, Parameter<double> > p_doubles_;
        std::unordered_map<std::string, Parameter<int> > p_ints_;
        std::unordered_map<std::string, Parameter<bool> > p_bools_;

        void init( Parameters* );
        std::unique_ptr<ParametersList> proc_params_;
        std::string proc_name_, hadr_name_, integr_type_;
    };

    //----- specialised registerers

    template<> inline void LpairHandler::registerParameter<std::string>( const char* key, const char* description, std::string* def ) { p_strings_.insert( std::make_pair( key, Parameter<std::string>( key, description, def ) ) ); }
    template<> inline void LpairHandler::registerParameter<double>( const char* key, const char* description, double* def ) { p_doubles_.insert( std::make_pair( key, Parameter<double>( key, description, def ) ) ); }
    template<> inline void LpairHandler::registerParameter<int>( const char* key, const char* description, int* def ) { p_ints_.insert( std::make_pair( key, Parameter<int>( key, description, def ) ) ); }
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

