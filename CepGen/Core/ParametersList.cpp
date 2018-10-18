#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  ParametersList&
  ParametersList::operator+=( const ParametersList& oth )
  {
    param_values_.insert( oth.param_values_.begin(), oth.param_values_.end() );
    int_values_.insert( oth.int_values_.begin(), oth.int_values_.end() );
    dbl_values_.insert( oth.dbl_values_.begin(), oth.dbl_values_.end() );
    str_values_.insert( oth.str_values_.begin(), oth.str_values_.end() );
    vec_param_values_.insert( oth.vec_param_values_.begin(), oth.vec_param_values_.end() );
    vec_int_values_.insert( oth.vec_int_values_.begin(), oth.vec_int_values_.end() );
    vec_dbl_values_.insert( oth.vec_dbl_values_.begin(), oth.vec_dbl_values_.end() );
    vec_str_values_.insert( oth.vec_str_values_.begin(), oth.vec_str_values_.end() );
    return *this;
  }

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

  template<typename T> T&
  ParametersList::operator[]( std::string key )
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

  template<> ParametersList&
  ParametersList::operator[]<ParametersList>( std::string key )
  {
    for ( auto& kv : param_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return param_values_[key];
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

  template<> std::vector<ParametersList>&
  ParametersList::operator[]<std::vector<ParametersList> >( std::string key )
  {
    for ( auto& kv : vec_param_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return vec_param_values_[key];
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

  template<> int&
  ParametersList::operator[]<int>( std::string key )
  {
    for ( auto& kv : int_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return int_values_[key];
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

  template<> std::vector<int>&
  ParametersList::operator[]<std::vector<int> >( std::string key )
  {
    for ( auto& kv : vec_int_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return vec_int_values_[key];
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

  template<> double&
  ParametersList::operator[]<double>( std::string key )
  {
    for ( auto& kv : dbl_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return dbl_values_[key];
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

  template<> std::vector<double>&
  ParametersList::operator[]<std::vector<double> >( std::string key )
  {
    for ( auto& kv : vec_dbl_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return vec_dbl_values_[key];
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

  template<> std::string&
  ParametersList::operator[]<std::string>( std::string key )
  {
    for ( auto& kv : str_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return str_values_[key];
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

  template<> std::vector<std::string>&
  ParametersList::operator[]<std::vector<std::string> >( std::string key )
  {
    for ( auto& kv : vec_str_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return vec_str_values_[key];
  }

  template<> void
  ParametersList::set<std::vector<std::string> >( std::string key, const std::vector<std::string>& value )
  {
    vec_str_values_[key] = value;
  }

  //------------------------------------------------------------------
  // limits-type attributes
  //------------------------------------------------------------------

  template<> Limits
  ParametersList::get<Limits>( std::string key, Limits def ) const
  {
    for ( const auto& kv : lim_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> Limits&
  ParametersList::operator[]<Limits>( std::string key )
  {
    for ( auto& kv : lim_values_ )
      if ( kv.first.compare( key ) == 0 )
        return kv.second;
    return lim_values_[key];
  }

  template<> void
  ParametersList::set<Limits>( std::string key, const Limits& value )
  {
    lim_values_[key] = value;
  }
}

