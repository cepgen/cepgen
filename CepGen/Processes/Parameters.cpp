#include "CepGen/Processes/Parameters.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace Process
  {
    std::ostream&
    operator<<( std::ostream& os, const Parameters& params )
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
    Parameters::get( std::string key, T def ) const
    {
      throw CG_FATAL( "ProcessParameters" ) << "Invalid type retrieved for key=" << key << "!";
    }

    template<typename T> void
    Parameters::set( std::string key, const T& value )
    {
      throw CG_FATAL( "ProcessParameters" ) << "Invalid type to be set for key=" << key << "!";
    }

    //------------------------------------------------------------------
    // sub-parameters-type attributes
    //------------------------------------------------------------------

    template<> Parameters
    Parameters::get<Parameters>( std::string key, Parameters def ) const
    {
      for ( const auto& kv : param_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<Parameters>( std::string key, const Parameters& value )
    {
      param_values_[key] = value;
    }

    template<> std::vector<Parameters>
    Parameters::get<std::vector<Parameters> >( std::string key, std::vector<Parameters> def ) const
    {
      for ( const auto& kv : vec_param_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<std::vector<Parameters> >( std::string key, const std::vector<Parameters>& value )
    {
      vec_param_values_[key] = value;
    }

    //------------------------------------------------------------------
    // integer-type attributes
    //------------------------------------------------------------------

    template<> int
    Parameters::get<int>( std::string key, int def ) const
    {
      for ( const auto& kv : int_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<int>( std::string key, const int& value )
    {
      int_values_[key] = value;
    }

    template<> std::vector<int>
    Parameters::get<std::vector<int> >( std::string key, std::vector<int> def ) const
    {
      for ( const auto& kv : vec_int_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<std::vector<int> >( std::string key, const std::vector<int>& value )
    {
      vec_int_values_[key] = value;
    }

    //------------------------------------------------------------------
    // floating point-type attributes
    //------------------------------------------------------------------

    template<> double
    Parameters::get<double>( std::string key, double def ) const
    {
      for ( const auto& kv : dbl_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<double>( std::string key, const double& value )
    {
      dbl_values_[key] = value;
    }

    template<> std::vector<double>
    Parameters::get<std::vector<double> >( std::string key, std::vector<double> def ) const
    {
      for ( const auto& kv : vec_dbl_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<std::vector<double> >( std::string key, const std::vector<double>& value )
    {
      vec_dbl_values_[key] = value;
    }

    //------------------------------------------------------------------
    // string-type attributes
    //------------------------------------------------------------------

    template<> std::string
    Parameters::get<std::string>( std::string key, std::string def ) const
    {
      for ( const auto& kv : str_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<std::string>( std::string key, const std::string& value )
    {
      str_values_[key] = value;
    }

    template<> std::vector<std::string>
    Parameters::get<std::vector<std::string> >( std::string key, std::vector<std::string> def ) const
    {
      for ( const auto& kv : vec_str_values_ )
        if ( kv.first.compare( key ) == 0 )
          return kv.second;
      CG_DEBUG( "ProcessParameters" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }

    template<> void
    Parameters::set<std::vector<std::string> >( std::string key, const std::vector<std::string>& value )
    {
      vec_str_values_[key] = value;
    }
  }
}
