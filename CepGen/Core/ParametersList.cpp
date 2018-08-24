#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  std::ostream&
  operator<<( std::ostream& os, const ParametersList& params )
  {
    for ( const auto& kv : params.int_values_ )
      os << "\n" << kv.first << ": int(" << kv.second << ")";
    for ( const auto& kv : params.dbl_values_ )
      os << "\n" << kv.first << ": double(" << kv.second << ")";
    for ( const auto& kv : params.str_values_ )
      os << "\n" << kv.first << ": string(" << kv.second << ")";
    for ( const auto& kv : params.param_values_ )
      os << "\n" << kv.first << ": param({" << kv.second << "\n})";
    for ( const auto& kv : params.vec_int_values_ ) {
      os << "\n" << kv.first << ": vint(";
      bool first = true;
      for ( const auto& v : kv.second ) {
        os << ( first ? "" : ", " ) << v;
        first = false;
      }
      os << ")";
    }
    for ( const auto& kv : params.vec_dbl_values_ ) {
      os << "\n" << kv.first << ": vdouble(";
      bool first = true;
      for ( const auto& v : kv.second ) {
        os << ( first ? "" : ", " ) << v;
        first = false;
      }
      os << ")";
    }
    for ( const auto& kv : params.vec_str_values_ ) {
      os << "\n" << kv.first << ": vstring(";
      bool first = true;
      for ( const auto& v : kv.second ) {
        os << ( first ? "" : ", " ) << v;
        first = false;
      }
      os << ")";
    }
    return os;
  }

  //------------------------------------------------------------------
  // default template (placeholders)
  //------------------------------------------------------------------

  template<typename T> T
  ParametersList::get( std::string key, T def ) const
  {
    throw CG_FATAL( "ParametersList" ) << "Invalid type retrieved for key=" << key << "!";
  }

  template<typename T> void
  ParametersList::set( std::string key, const T& value )
  {
    throw CG_FATAL( "ParametersList" ) << "Invalid type to be set for key=" << key << "!";
  }

  //------------------------------------------------------------------
  // sub-parameters-type attributes
  //------------------------------------------------------------------

  template<> ParametersList
  ParametersList::get<ParametersList>( std::string key, ParametersList def ) const
  {
    for ( const auto& kv : param_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<ParametersList>( std::string key, const ParametersList& value )
  {
    param_values_[key] = value;
  }

  template<> std::vector<ParametersList>
  ParametersList::get<std::vector<ParametersList> >( std::string key, std::vector<ParametersList> def ) const
  {
    for ( const auto& kv : vec_param_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<std::vector<ParametersList> >( std::string key, const std::vector<ParametersList>& value )
  {
    vec_param_values_[key] = value;
  }

  //------------------------------------------------------------------
  // integer-type attributes
  //------------------------------------------------------------------

  template<> int
  ParametersList::get<int>( std::string key, int def ) const
  {
    for ( const auto& kv : int_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<int>( std::string key, const int& value )
  {
    int_values_[key] = value;
  }

  template<> std::vector<int>
  ParametersList::get<std::vector<int> >( std::string key, std::vector<int> def ) const
  {
    for ( const auto& kv : vec_int_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<std::vector<int> >( std::string key, const std::vector<int>& value )
  {
    vec_int_values_[key] = value;
  }

  //------------------------------------------------------------------
  // floating point-type attributes
  //------------------------------------------------------------------

  template<> double
  ParametersList::get<double>( std::string key, double def ) const
  {
    for ( const auto& kv : dbl_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<double>( std::string key, const double& value )
  {
    dbl_values_[key] = value;
  }

  template<> std::vector<double>
  ParametersList::get<std::vector<double> >( std::string key, std::vector<double> def ) const
  {
    for ( const auto& kv : vec_dbl_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<std::vector<double> >( std::string key, const std::vector<double>& value )
  {
    vec_dbl_values_[key] = value;
  }

  //------------------------------------------------------------------
  // string-type attributes
  //------------------------------------------------------------------

  template<> std::string
  ParametersList::get<std::string>( std::string key, std::string def ) const
  {
    for ( const auto& kv : str_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<std::string>( std::string key, const std::string& value )
  {
    str_values_[key] = value;
  }

  template<> std::vector<std::string>
  ParametersList::get<std::vector<std::string> >( std::string key, std::vector<std::string> def ) const
  {
    for ( const auto& kv : vec_str_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> void
  ParametersList::set<std::vector<std::string> >( std::string key, const std::vector<std::string>& value )
  {
    vec_str_values_[key] = value;
  }
}
