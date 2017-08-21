#ifndef CepGen_Cards_LpairReader_h
#define CepGen_Cards_LpairReader_h

#include "Handler.h"

#include <fstream>
#include <string>
#include <map>

namespace CepGen
{
  namespace Cards
  {
    class LpairReader : public Handler
    {
      public:
        LpairReader( const char* file );

        void store( const char* file );

      private:
        template<class T> struct Parameter {
          Parameter( const char* key, const char* descr, T* value ) : key( key ), description( descr ), value( value ) {}
          std::string key, description;
          T* value;
        };
        template<class T> void registerParameter( const char* key, const char* description, T* def ) {}
        template<class T> void setValue( const char* key, T value ) {}
        template<class T> T getValue( const char* key ) const {}

        void setParameter( std::string key, std::string value );
        std::string getParameter( std::string key ) const;
        std::string getDescription( std::string key ) const;

        std::map<std::string, Parameter<std::string> > p_strings_;
        std::map<std::string, Parameter<double> > p_doubles_;
        std::map<std::string, Parameter<unsigned int> > p_ints_;
        std::map<std::string, Parameter<bool> > p_bools_;

        void init( Parameters* );
        std::string proc_name_, hadr_name_;
    };

    //----- specialised registerers

    template<> inline void LpairReader::registerParameter<std::string>( const char* key, const char* description, std::string* def ) { p_strings_.insert( std::make_pair( key, Parameter<std::string>( key, description, def ) ) ); }
    template<> inline void LpairReader::registerParameter<double>( const char* key, const char* description, double* def ) { p_doubles_.insert( std::make_pair( key, Parameter<double>( key, description, def ) ) ); }
    template<> inline void LpairReader::registerParameter<unsigned int>( const char* key, const char* description, unsigned int* def ) { p_ints_.insert( std::make_pair( key, Parameter<unsigned int>( key, description, def ) ) ); }
    template<> inline void LpairReader::registerParameter<bool>( const char* key, const char* description, bool* def ) { p_bools_.insert( std::make_pair( key, Parameter<bool>( key, description, def ) ) ); }

    //----- specialised setters

    template<> inline void LpairReader::setValue<std::string>( const char* key, std::string value ) {
      auto it = p_strings_.find( key );
      if ( it != p_strings_.end() ) *it->second.value = value;
    }
    template<> inline void LpairReader::setValue<double>( const char* key, double value ) {
      auto it = p_doubles_.find( key );
      if ( it != p_doubles_.end() ) *it->second.value = value;
    }
    template<> inline void LpairReader::setValue<unsigned int>( const char* key, unsigned int value ) {
      auto it = p_ints_.find( key );
      if ( it != p_ints_.end() ) *it->second.value = value;
    }
    template<> inline void LpairReader::setValue<bool>( const char* key, bool value ) {
      auto it = p_bools_.find( key );
      if ( it != p_bools_.end() ) *it->second.value = value;
    }

    //----- specialised getters

    template<> inline std::string LpairReader::getValue( const char* key ) const {
      const auto& it = p_strings_.find( key );
      if ( it != p_strings_.end() ) return *it->second.value;
      return "null";
    }
    template<> inline double LpairReader::getValue( const char* key ) const {
      const auto& it = p_doubles_.find( key );
      if ( it != p_doubles_.end() ) return *it->second.value;
      return -999.;
    }
    template<> inline unsigned int LpairReader::getValue( const char* key ) const {
      const auto& it = p_ints_.find( key );
      if ( it != p_ints_.end() ) return *it->second.value;
      return 999;
    }
    template<> inline bool LpairReader::getValue( const char* key ) const {
      const auto& it = p_bools_.find( key );
      if ( it != p_bools_.end() ) return *it->second.value;
      return true;
    }
  }
}

#endif
