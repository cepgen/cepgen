#ifndef CepGen_Core_ParametersList_h
#define CepGen_Core_ParametersList_h

#include "CepGen/Physics/Limits.h"
#include "CepGen/Core/Hasher.h"

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

namespace cepgen
{
  /// Parameters container
  class ParametersList
  {
    private:
      template<typename T> struct default_arg
      {
        static T get() { return T(); }
      };

    public:
      ParametersList() = default;
      ~ParametersList() = default; // required for unique_ptr initialisation!
      /// Check if a given parameter is handled in this list
      template<typename T> bool has( std::string key ) const;
      /// Get a parameter value
      template<typename T> T get( std::string key, const T& def = default_arg<T>::get() ) const;
      /// Reference to a parameter value
      template<typename T> T& operator[]( std::string key );
      /// Set a parameter value
      template<typename T> ParametersList& set( std::string key, const T& value );
      /// Concatenate two parameters containers
      ParametersList& operator+=( const ParametersList& oth );

      /// List of keys handled in this list of parameters
      std::vector<std::string> keys() const;
      /// Get a string-converted version of a value
      std::string getString( const std::string& key ) const;

      /// Human-readable version of a parameters container
      friend std::ostream& operator<<( std::ostream& os, const ParametersList& );

    private:
      std::map<std::string,ParametersList> param_values_;
      std::unordered_map<std::string,int> int_values_;
      std::unordered_map<std::string,double> dbl_values_;
      std::unordered_map<std::string,std::string> str_values_;
      std::unordered_map<std::string,Limits> lim_values_;
      std::unordered_map<std::string,std::vector<ParametersList> > vec_param_values_;
      std::unordered_map<std::string,std::vector<int> > vec_int_values_;
      std::unordered_map<std::string,std::vector<double> > vec_dbl_values_;
      std::unordered_map<std::string,std::vector<std::string> > vec_str_values_;
  };
  /// Check if an integer parameter is handled
  template<> inline bool ParametersList::has<int>( std::string key ) const { return int_values_.count( key ) != 0; }
  /// Get an integer parameter value
  template<> int ParametersList::get<int>( std::string key, const int& def ) const;
  /// Reference to an integer parameter value
  template<> int& ParametersList::operator[]<int>( std::string key );
  /// Set an integer parameter value
  template<> ParametersList& ParametersList::set<int>( std::string key, const int& value );
  /// Check if a vector of integers parameter is handled
  template<> inline bool ParametersList::has<std::vector<int> >( std::string key ) const { return vec_int_values_.count( key ) != 0; }
  /// Get a vector of integers parameter value
  template<> std::vector<int> ParametersList::get<std::vector<int> >( std::string key, const std::vector<int>& def ) const;
  /// Reference to a vector of integers parameter value
  template<> std::vector<int>& ParametersList::operator[]<std::vector<int> >( std::string key );
  /// Set a vector of integers parameter value
  template<> ParametersList& ParametersList::set<std::vector<int> >( std::string key, const std::vector<int>& value );

  /// Check if a boolean parameter is handled
  template<> inline bool ParametersList::has<bool>( std::string key ) const { return has<int>( key ); }
  /// Get a boolean parameter value
  template<> inline bool ParametersList::get<bool>( std::string key, const bool& def ) const { return static_cast<bool>( get<int>( key, def ) ); }
  /// Reference to a boolean parameter value
  template<> inline bool& ParametersList::operator[]<bool>( std::string key ) { return (bool&)operator[]<int>( key ); }
  /// Set a boolean parameter value
  template<> inline ParametersList& ParametersList::set<bool>( std::string key, const bool& value ) { return set<int>( key, static_cast<bool>( value ) ); }

  /// Check if a double floating point parameter is handled
  template<> inline bool ParametersList::has<double>( std::string key ) const { return dbl_values_.count( key ) != 0; }
  /// Get a double floating point parameter value
  template<> double ParametersList::get<double>( std::string key, const double& def ) const;
  /// Reference to a double floating point parameter value
  template<> double& ParametersList::operator[]<double>( std::string key );
  /// Set a double floating point parameter value
  template<> ParametersList& ParametersList::set<double>( std::string key, const double& value );
  /// Check if a vector of double floating point parameter is handled
  template<> inline bool ParametersList::has<std::vector<double> >( std::string key ) const { return vec_dbl_values_.count( key ) != 0; }
  /// Get a vector of double floating point parameter value
  template<> std::vector<double> ParametersList::get<std::vector<double> >( std::string key, const std::vector<double>& def ) const;
  /// Reference to a vector of double floating point parameter value
  template<> std::vector<double>& ParametersList::operator[]<std::vector<double> >( std::string key );
  /// Set a vector of double floating point parameter value
  template<> ParametersList& ParametersList::set<std::vector<double> >( std::string key, const std::vector<double>& value );

  /// Check if a string parameter is handled
  template<> inline bool ParametersList::has<std::string>( std::string key ) const { return str_values_.count( key ) != 0; }
  /// Get a string parameter value
  template<> std::string ParametersList::get<std::string>( std::string key, const std::string& def ) const;
  /// Reference to a string parameter value
  template<> std::string& ParametersList::operator[]<std::string>( std::string key );
  /// Set a string parameter value
  template<> ParametersList& ParametersList::set<std::string>( std::string key, const std::string& value );
  /// Check if a vector of strings parameter is handled
  template<> inline bool ParametersList::has<std::vector<std::string> >( std::string key ) const { return vec_str_values_.count( key ) != 0; }
  /// Get a vector of strings parameter value
  template<> std::vector<std::string> ParametersList::get<std::vector<std::string> >( std::string key, const std::vector<std::string>& def ) const;
  /// Reference to a vector of strings parameter value
  template<> std::vector<std::string>& ParametersList::operator[]<std::vector<std::string> >( std::string key );
  /// Set a vector of strings parameter value
  template<> ParametersList& ParametersList::set<std::vector<std::string> >( std::string key, const std::vector<std::string>& value );

  /// Check if a limits parameter is handled
  template<> inline bool ParametersList::has<Limits>( std::string key ) const { return lim_values_.count( key ) != 0; }
  /// Get a boundary limits parameter value
  template<> Limits ParametersList::get<Limits>( std::string key, const Limits& def ) const;
  /// Reference to a boundary limits parameter value
  template<> Limits& ParametersList::operator[]<Limits>( std::string key );
  /// Set a boundary limits parameter value
  template<> ParametersList& ParametersList::set<Limits>( std::string key, const Limits& value );

  /// Check if a parameters list parameter is handled
  template<> inline bool ParametersList::has<ParametersList>( std::string key ) const { return param_values_.count( key ) != 0; }
  /// Get a parameters list parameter value
  template<> ParametersList ParametersList::get<ParametersList>( std::string key, const ParametersList& def ) const;
  /// Reference to a parameters list parameter value
  template<> ParametersList& ParametersList::operator[]<ParametersList>( std::string key );
  /// Set a parameters list parameter value
  template<> ParametersList& ParametersList::set<ParametersList>( std::string key, const ParametersList& value );
  /// Check if a vector of parameters lists is handled
  template<> inline bool ParametersList::has<std::vector<ParametersList> >( std::string key ) const { return vec_param_values_.count( key ) != 0; }
  /// Get a vector of parameters list parameter value
  template<> std::vector<ParametersList> ParametersList::get<std::vector<ParametersList> >( std::string key, const std::vector<ParametersList>& def ) const;
  /// Reference to a vector of parameters list parameter value
  template<> std::vector<ParametersList>& ParametersList::operator[]<std::vector<ParametersList> >( std::string key );
  /// Set a vector of parameters list parameter value
  template<> ParametersList& ParametersList::set<std::vector<ParametersList> >( std::string key, const std::vector<ParametersList>& value );
}

#endif
