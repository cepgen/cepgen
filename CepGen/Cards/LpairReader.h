#ifndef CepGen_Cards_LpairReader_h
#define CepGen_Cards_LpairReader_h

#include "Handler.h"

#include <fstream>
#include <string>

namespace CepGen
{
  namespace Cards
  {
    class LpairReader : public Handler
    {
      public:
        LpairReader( const char* file );

        void store( const char* file ) const;

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

        std::vector<Parameter<std::string> > p_strings_;
        std::vector<Parameter<double> > p_doubles_;
        std::vector<Parameter<unsigned int> > p_ints_;
        std::vector<Parameter<bool> > p_bools_;

        //void parseKinematics( const libconfig::Setting& );
        std::map<std::string,std::string> params_map_;
    };
    template<> inline void LpairReader::registerParameter<std::string>( const char* key, const char* description, std::string* def ) { p_strings_.emplace_back( key, description, def ); }
    template<> inline void LpairReader::registerParameter<double>( const char* key, const char* description, double* def ) { p_doubles_.emplace_back( key, description, def ); }
    template<> inline void LpairReader::registerParameter<unsigned int>( const char* key, const char* description, unsigned int* def ) { p_ints_.emplace_back( key, description, def ); }
    template<> inline void LpairReader::registerParameter<bool>( const char* key, const char* description, bool* def ) { p_bools_.emplace_back( key, description, def ); }
    template<> inline void LpairReader::setValue<std::string>( const char* key, std::string value ) {
      for ( std::vector<Parameter<std::string> >::iterator it = p_strings_.begin(); it != p_strings_.end(); ++it ) {
        if ( it->key == key ) *it->value = value;
      }
    }
    template<> inline void LpairReader::setValue<double>( const char* key, double value ) {
      for ( std::vector<Parameter<double> >::iterator it = p_doubles_.begin(); it != p_doubles_.end(); ++it ) {
        if ( it->key == key ) *it->value = value;
      }
    }
    template<> inline void LpairReader::setValue<unsigned int>( const char* key, unsigned int value ) {
      for ( std::vector<Parameter<unsigned int> >::iterator it = p_ints_.begin(); it != p_ints_.end(); ++it ) {
        if ( it->key == key ) *it->value = value;
      }
    }
    template<> inline void LpairReader::setValue<bool>( const char* key, bool value ) {
      for ( std::vector<Parameter<bool> >::iterator it = p_bools_.begin(); it != p_bools_.end(); ++it ) {
        if ( it->key == key ) *it->value = value;
      }
    }
    template<> inline std::string LpairReader::getValue( const char* key ) const {
      for ( std::vector<Parameter<std::string> >::const_iterator it = p_strings_.begin(); it != p_strings_.end(); ++it ) {
        if ( it->key == key ) return *it->value;
      }
      return "null";
    }
    template<> inline double LpairReader::getValue( const char* key ) const {
      for ( std::vector<Parameter<double> >::const_iterator it = p_doubles_.begin(); it != p_doubles_.end(); ++it ) {
        if ( it->key == key ) return *it->value;
      }
      return -999.;
    }
    template<> inline unsigned int LpairReader::getValue( const char* key ) const {
      for ( std::vector<Parameter<unsigned int> >::const_iterator it = p_ints_.begin(); it != p_ints_.end(); ++it ) {
        if ( it->key == key ) return *it->value;
      }
      return 999;
    }
  }
}

#endif
