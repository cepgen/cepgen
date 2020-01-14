#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Physics/PDG.h"

#include <iomanip>

namespace cepgen
{
  ParametersList::ParametersList( const ParametersList& oth ) :
    param_values_( oth.param_values_ ), int_values_( oth.int_values_ ),
    dbl_values_( oth.dbl_values_ ), str_values_( oth.str_values_ ),
    lim_values_( oth.lim_values_ ),
    vec_param_values_( oth.vec_param_values_ ), vec_int_values_( oth.vec_int_values_ ),
    vec_dbl_values_( oth.vec_dbl_values_ ), vec_str_values_( oth.vec_str_values_ )
  {}

  ParametersList&
  ParametersList::operator+=( const ParametersList& oth )
  {
    param_values_.insert( oth.param_values_.begin(), oth.param_values_.end() );
    int_values_.insert( oth.int_values_.begin(), oth.int_values_.end() );
    dbl_values_.insert( oth.dbl_values_.begin(), oth.dbl_values_.end() );
    str_values_.insert( oth.str_values_.begin(), oth.str_values_.end() );
    lim_values_.insert( oth.lim_values_.begin(), oth.lim_values_.end() );
    vec_param_values_.insert( oth.vec_param_values_.begin(), oth.vec_param_values_.end() );
    vec_int_values_.insert( oth.vec_int_values_.begin(), oth.vec_int_values_.end() );
    vec_dbl_values_.insert( oth.vec_dbl_values_.begin(), oth.vec_dbl_values_.end() );
    vec_str_values_.insert( oth.vec_str_values_.begin(), oth.vec_str_values_.end() );
    return *this;
  }

  bool
  ParametersList::empty() const
  {
    return keys().empty()
      || keys() == std::vector<std::string>{ MODULE_NAME };
  }

  std::ostream&
  operator<<( std::ostream& os, const ParametersList& params )
  {
    const auto beg = os.tellp();
    for ( const auto& kv : params.int_values_ )
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=int(" << kv.second << ")";
    for ( const auto& kv : params.dbl_values_ )
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=double(" << kv.second << ")";
    for ( const auto& kv : params.str_values_ )
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=string(" << kv.second << ")";
    for ( const auto& kv : params.param_values_ )
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=param({" << kv.second << "})";
    for ( const auto& kv : params.lim_values_ )
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=limits(" << kv.second << ")";
    for ( const auto& kv : params.vec_int_values_ ) {
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=vint(";
      bool first = true;
      for ( const auto& v : kv.second )
        os << ( first ? "" : ", " ) << v, first = false;
      os << ")";
    }
    for ( const auto& kv : params.vec_dbl_values_ ) {
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=vdouble(";
      bool first = true;
      for ( const auto& v : kv.second )
        os << ( first ? "" : ", " ) << v, first = false;
      os << ")";
    }
    for ( const auto& kv : params.vec_str_values_ ) {
      os << ( os.tellp() > beg ? ", " : "" ) << kv.first << "=vstring(";
      bool first = true;
      for ( const auto& v : kv.second )
        os << ( first ? "" : ", " ) << v, first = false;
      os << ")";
    }
    return os;
  }

  std::vector<std::string>
  ParametersList::keys() const
  {
    std::vector<std::string> out;
    for ( const auto& p :     param_values_ ) out.emplace_back( p.first );
    for ( const auto& p : vec_param_values_ ) out.emplace_back( p.first );
    for ( const auto& p :       int_values_ ) out.emplace_back( p.first );
    for ( const auto& p :   vec_int_values_ ) out.emplace_back( p.first );
    for ( const auto& p :       dbl_values_ ) out.emplace_back( p.first );
    for ( const auto& p :   vec_dbl_values_ ) out.emplace_back( p.first );
    for ( const auto& p :       str_values_ ) out.emplace_back( p.first );
    for ( const auto& p :   vec_str_values_ ) out.emplace_back( p.first );
    for ( const auto& p :       lim_values_ ) out.emplace_back( p.first );
    return out;
  }

  std::string
  ParametersList::getString( const std::string& key ) const
  {
    std::ostringstream os;
    if ( has<ParametersList>( key ) )   os << "params{" << get<ParametersList>( key ) << "}";
    else if ( has<int>( key ) )         os << get<int>( key );
    else if ( has<double>( key ) )      os << get<double>( key );
    else if ( has<std::string>( key ) ) os << get<std::string>( key );
    else if ( has<Limits>( key ) )      os << get<Limits>( key );
    else if ( has<std::vector<ParametersList> >( key ) ) {
      bool first = true;
      for ( const auto& p : get<std::vector<ParametersList> >( key ) ) {
        os << ( first ? "" : ", " ) << p;
        first = false;
      }
    }
    else if ( has<std::vector<int> >( key ) ) {
      bool first = true;
      for ( const auto& p : get<std::vector<int> >( key ) ) {
        os << ( first ? "" : ", " ) << p;
        first = false;
      }
    }
    else if ( has<std::vector<double> >( key ) ) {
      bool first = true;
      for ( const auto& p : get<std::vector<double> >( key ) ) {
        os << ( first ? "" : ", " ) << p;
        first = false;
      }
    }
    else if ( has<std::vector<std::string> >( key ) ) {
      bool first = true;
      for ( const auto& p : get<std::vector<std::string> >( key ) ) {
        os << ( first ? "" : ", " ) << p;
        first = false;
      }
    }
    return os.str();
  }

  //------------------------------------------------------------------
  // default template (placeholders)
  //------------------------------------------------------------------

  template<typename T> bool
  ParametersList::has( std::string key ) const
  {
    throw CG_FATAL( "ParametersList" ) << "Invalid type for key=" << key << "!";
  }

  template<typename T> T
  ParametersList::get( std::string key, const T& def ) const
  {
    throw CG_FATAL( "ParametersList" ) << "Invalid type retrieved for key=" << key << "!";
  }

  template<typename T> T&
  ParametersList::operator[]( std::string key )
  {
    throw CG_FATAL( "ParametersList" ) << "Invalid type retrieved for key=" << key << "!";
  }

  template<typename T> ParametersList&
  ParametersList::set( std::string key, const T& value )
  {
    throw CG_FATAL( "ParametersList" ) << "Invalid type to be set for key=" << key << "!";
  }

  //------------------------------------------------------------------
  // sub-parameters-type attributes
  //------------------------------------------------------------------

  template<> ParametersList
  ParametersList::get<ParametersList>( std::string key, const ParametersList& def ) const
  {
    for ( const auto& kv : param_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> std::vector<ParametersList>
  ParametersList::get<std::vector<ParametersList> >( std::string key, const std::vector<ParametersList>& def ) const
  {
    for ( const auto& kv : vec_param_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  //------------------------------------------------------------------
  // integer-type attributes
  //------------------------------------------------------------------

  template<> int
  ParametersList::get<int>( std::string key, const int& def ) const
  {
    for ( const auto& kv : int_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> std::vector<int>
  ParametersList::get<std::vector<int> >( std::string key, const std::vector<int>& def ) const
  {
    for ( const auto& kv : vec_int_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  //------------------------------------------------------------------
  // floating point-type attributes
  //------------------------------------------------------------------

  template<> double
  ParametersList::get<double>( std::string key, const double& def ) const
  {
    for ( const auto& kv : dbl_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> std::vector<double>
  ParametersList::get<std::vector<double> >( std::string key, const std::vector<double>& def ) const
  {
    for ( const auto& kv : vec_dbl_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  //------------------------------------------------------------------
  // string-type attributes
  //------------------------------------------------------------------

  template<> std::string
  ParametersList::get<std::string>( std::string key, const std::string& def ) const
  {
    for ( const auto& kv : str_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  template<> std::vector<std::string>
  ParametersList::get<std::vector<std::string> >( std::string key, const std::vector<std::string>& def ) const
  {
    for ( const auto& kv : vec_str_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  //------------------------------------------------------------------
  // limits-type attributes
  //------------------------------------------------------------------

  template<> Limits
  ParametersList::get<Limits>( std::string key, const Limits& def ) const
  {
    for ( const auto& kv : lim_values_ )
      if ( kv.first == key )
        return kv.second;
    CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
    return def;
  }

  //------------------------------------------------------------------
  // particle properties-type attributes
  //------------------------------------------------------------------

  template<> ParticleProperties
  ParametersList::get<ParticleProperties>( std::string key, const ParticleProperties& def ) const
  {
    if ( has<ParametersList>( key ) ) {
      const auto& plist = get<ParametersList>( key );
      ParticleProperties out;
      const auto& pdgid = plist.get<int>( "pdgid", 0 );
      if ( PDG::get().has( pdgid ) )
        out = PDG::get()( pdgid );
      else
        out.pdgid = pdgid;
      bool modified = false;
      if ( plist.has<std::string>( "name" ) )
        out.name = plist.get<std::string>( "name" ), modified = true;
      if ( plist.has<std::string>( "description" ) )
        out.description = plist.get<std::string>( "description" ), modified = true;
      if ( plist.has<int>( "colours" ) )
        out.colours = plist.get<int>( "colours" ), modified = true;
      if ( plist.has<double>( "mass" ) )
        out.mass = plist.get<double>( "mass" ), modified = true;
      if ( plist.has<double>( "width" ) )
        out.width = plist.get<double>( "width" ), modified = true;
      if ( plist.has<double>( "charge" ) )
        out.charge = short( plist.get<double>( "charge" )*3 ), modified = true;
      if ( plist.has<bool>( "fermion" ) )
        out.fermion = plist.get<bool>( "fermion" ), modified = true;
      if ( modified )
        PDG::get().define( out );
      return out;
    }
    else if ( has<int>( key ) )
      return PDG::get()( get<int>( key ) );
    else {
      CG_DEBUG( "ParametersList" ) << "Failed to retrieve parameter with key=" << key << ".";
      return def;
    }
  }

  template<> ParametersList&
  ParametersList::set<ParticleProperties>( std::string key, const ParticleProperties& value )
  {
    return set<ParametersList>( key, ParametersList()
      .set<int>( "pdgid", value.pdgid )
      .set<std::string>( "name", value.name )
      .set<std::string>( "description", value.description )
      .set<int>( "colours", value.colours )
      .set<double>( "mass", value.mass )
      .set<double>( "width", value.width )
      .set<double>( "charge", value.charge*1./3 )
      .set<bool>( "fermion", value.fermion ) );
  }
}
