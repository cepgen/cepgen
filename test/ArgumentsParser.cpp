#include "ArgumentsParser.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <sstream>
#include <algorithm>
#include <cstring>

namespace cepgen
{
  ArgumentsParser::ArgumentsParser( int argc, char* argv[] ) :
    help_str_( { { "help", 'h' } } )
  {
    command_name_ = argv[0];
    //--- first remove the program name
    std::vector<std::string> args_tmp;
    if ( argc > 1 ) {
      args_tmp.resize( argc-1 );
      std::copy( argv+1, argv+argc, args_tmp.begin() );
    }
    //--- then build the arguments list
    for ( const auto& arg : args_tmp ) {
      //--- skip the '='
      const auto eq_pos = arg.find( '=' );
      if ( eq_pos != std::string::npos ) {
        args_.emplace_back( arg.substr( 0, eq_pos ) );
        args_.emplace_back( arg.substr( eq_pos+1 ) );
        continue;
      }
      args_.emplace_back( arg );
    }
  }

  void
  ArgumentsParser::print_help() const
  {
    CG_INFO( "ArgumentsParser" ) << help_message();
    exit( 0 );
  }

  void
  ArgumentsParser::dump() const
  {
    std::ostringstream os;
    os
      << "List of parameters retrieved from command-line:";
    for ( const auto& par : params_ )
      os
        << "\n[" << par.name
        << ( par.sname != '\0' ? "|"+std::string( 1, par.sname ) : "" )
        << "] = " << par.value << ( par.optional ? ", optional" : "" );
    CG_INFO( "ArgumentsParser" ) << os.str();
  }

  ArgumentsParser&
  ArgumentsParser::parse()
  {
    if ( !args_.empty() )
      //--- check if help message is requested
      for ( const auto& str : help_str_ )
        if ( find( args_.begin(), args_.end(), "--"+str.name ) != args_.end()
          || find( args_.begin(), args_.end(), "-"+std::string( 1, str.sname ) ) != args_.end() )
          print_help();
    //--- loop over all parameters
    size_t i = 0;
    for ( auto& par : params_ ) {
      /*if ( !par.optional && args_.empty() )
        //--- if no arguments provided while at least one is required
        throw CG_FATAL( "ArgumentsParser" )
          << help_message()
          << " The following parameter was not set: '"
          << ( !par.name.empty() ? par.name : "<arg"+std::to_string( i )+">" )
          << "'.";*/
      if ( par.name.empty() ) {
        //--- no argument name ; fetching by index
        if ( i >= args_.size() )
          throw CG_FATAL( "ArgumentsParser" )
            << help_message()
            << " Failed to retrieve argument " << ( i+1 ) << " while required.";
        par.value = par.bool_variable ? "1" : args_.at( i );
      }
      else {
        const auto it_key = find( args_.begin(), args_.end(), "--"+par.name );
        const auto it_skey = find( args_.begin(), args_.end(), "-"+std::string( 1, par.sname ) );
        if ( it_key == args_.end() && it_skey == args_.end() ) { // not found
          if ( !par.optional )
            throw CG_FATAL( "ArgumentsParser" )
              << help_message()
              << " The following parameter was not set: '" << par.name << "'.";
        }
        const auto it_value = ( it_key != args_.end() )
          ? std::next( it_key )
          : std::next( it_skey );
        if ( it_value != args_.end() )
          par.value = *it_value;
        else if ( par.bool_variable )
          par.value = "1"; // if the flag is set, enabled by default
        else
          throw CG_FATAL( "ArgumentsParser" )
            << "Invalid value for parameter: " << par.name << ".";
      }
      par.parse();
      ++i;
    }
    return *this;
  }

  std::string
  ArgumentsParser::operator[]( std::string name ) const
  {
    for ( const auto& par : params_ ) {
      if ( "--"+par.name == name ) return par.value;
      if ( par.sname != '\0' && "-"+std::string( 1, par.sname ) == name )
        return par.value;
    }
    throw CG_FATAL( "ArgumentsParser" )
      << "The parameter \"" << name << "\" was not declared "
      << "in the arguments parser constructor!";
  }

  std::string
  ArgumentsParser::help_message() const
  {
    std::ostringstream oss;
    std::vector<std::pair<Parameter, size_t> > req_params, opt_params;
    oss << "Usage: " << command_name_;
    size_t i = 0;
    for ( const auto& par : params_ ) {
      if ( par.optional ) {
        opt_params.emplace_back( std::make_pair( par, i ) );
        oss << "[";
      }
      else
        req_params.emplace_back( std::make_pair( par, i ) );
      oss << ( !par.name.empty() ? " --" : " <arg"+std::to_string( i )+">" ) << par.name;
      if ( par.sname != '\0' ) oss << "|-" << par.sname;
      if ( par.optional )
        oss << "]";
      ++i;
    }
    if ( req_params.size() > 0 ) {
      oss << "\n    " << utils::s( "required argument", req_params.size() ) << ":";
      for ( const auto& par : req_params )
        oss << Form( ( par.first.sname != '\0' )
          ? "\n\t%s/-%1s\t%-28s"
          : "\n\t%s/%2s\t%-28s",
          ( !par.first.name.empty() ? "--"+par.first.name : "<arg"+std::to_string( par.second )+">" ).c_str(),
          &par.first.sname, par.first.description.c_str() );
    }
    if ( opt_params.size() > 0 ) {
      oss << "\n    " << utils::s( "optional argument", opt_params.size() ) << ":";
      for ( const auto& par : opt_params )
        oss << Form( ( par.first.sname != '\0' )
          ? "\n\t%s/-%1s\t%-28s\tdefault = '%s'"
          : "\n\t%s/%2s\t%-28s\tdefault = '%s'",
          ( !par.first.name.empty() ? "--"+par.first.name : "<arg"+std::to_string( par.second )+">" ).c_str(),
          &par.first.sname, par.first.description.c_str(),
          par.first.value.c_str() );
    }
    oss << std::endl;
    return oss.str();
  }

  //----- simple parameters

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::string default_value, std::string* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( default_value ), optional( true ),
    str_variable( var ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::string* var, char sname ) :
    Parameter( name, description, "", var, sname )
  {
    optional = false;
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, unsigned int default_value, unsigned int* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( std::to_string( default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( var ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, unsigned int* var, char sname ) :
    Parameter( name, description, 0, var, sname )
  {
    optional = false;
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, int default_value, int* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( Form( "%+i", default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( var ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, int* var, char sname ) :
    Parameter( name, description, 0, var, sname )
  {
    optional = false;
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, bool default_value, bool* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( Form( "%d", default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( var ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, bool* var, char sname ) :
    Parameter( name, description, 0, var, sname )
  {
    optional = false;
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, double default_value, double* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( Form( "%g", default_value ) ), optional( true ),
    str_variable( nullptr ), float_variable( var ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, double* var, char sname ) :
    Parameter( name, description, 0., var, sname )
  {
    optional = false;
  }

  //----- vector of parameters

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<std::string> default_value, std::vector<std::string>* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( "" ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( var ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {
    unsigned short i = 0;
    for ( const auto& str : default_value )
      value += ( ( ( i++ > 0 ) ? "," : "" )+str );
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<std::string>* var, char sname ) :
    Parameter( name, description, std::vector<std::string>{ { } }, var, sname )
  {
    optional = false;
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<int> default_value, std::vector<int>* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( "" ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( var ), vec_float_variable( nullptr )
  {
    unsigned short i = 0;
    for ( const auto& var : default_value )
      value += ( ( ( i++ > 0 ) ? "," : "" )+Form( "%d", var ) );
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<int>* var, char sname ) :
    Parameter( name, description, std::vector<int>{ { } }, var, sname )
  {
    optional = false;
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<double> default_value, std::vector<double>* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( "" ), optional( true ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( var )
  {
    unsigned short i = 0;
    for ( const auto& flt : default_value )
      value += ( ( ( i++ > 0 ) ? "," : "" )+Form( "%g", flt ) );
  }

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<double>* var, char sname ) :
    Parameter( name, description, std::vector<double>{ { } }, var, sname )
  {
    optional = false;
  }

  void
  ArgumentsParser::Parameter::parse()
  {
    if ( str_variable != nullptr )
      *str_variable = value;
    else if ( float_variable != nullptr )
      *float_variable = std::stod( value );
    else if ( int_variable != nullptr )
      *int_variable = std::stoi( value );
    else if ( uint_variable != nullptr )
      *uint_variable = std::stoi( value );
    else if ( bool_variable != nullptr ) {
      try {
        *bool_variable = ( std::stoi( value ) != 0 );
      } catch ( const std::invalid_argument& ) {
        *bool_variable = ( strcasecmp( "true", value.c_str() ) == 0
                        || strcasecmp( "yes", value.c_str() ) == 0 )
                        && strcasecmp( "false", value.c_str() ) != 0
                        && strcasecmp( "no", value.c_str() ) != 0;
      }
    }
    else if ( vec_str_variable != nullptr ) {
      std::istringstream iss( value ); std::string token;
      std::vector<std::string> vec_var;
      while ( std::getline( iss, token, ',' ) )
        vec_var.emplace_back( token );
      *vec_str_variable = vec_var;
    }
    else if ( vec_int_variable != nullptr ) {
      std::istringstream iss( value ); std::string token;
      std::vector<int> vec_var;
      while ( std::getline( iss, token, ',' ) )
        vec_var.emplace_back( std::stoi( token ) );
      *vec_int_variable = vec_var;
    }
    else if ( vec_float_variable != nullptr ) {
      std::istringstream iss( value ); std::string token;
      std::vector<double> vec_var;
      while ( std::getline( iss, token, ',' ) )
        vec_var.emplace_back( std::stod( token ) );
      *vec_float_variable = vec_var;
    }
  }
}

