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
      const auto eq_pos = arg.find( '=' );
      if ( eq_pos == std::string::npos ) {
        args_.emplace_back( arg );
        continue;
      }
      //--- skip the '='
      args_.emplace_back( arg.substr( 0, eq_pos ) );
      args_.emplace_back( arg.substr( eq_pos+1 ) );
    }
  }

  void
  ArgumentsParser::parse()
  {
    //--- check if help message is requested
    for ( const auto& str : help_str_ ) {
      if ( find( args_.begin(), args_.end(), "--"+str.name ) != args_.end()
        || find( args_.begin(), args_.end(), "-"+std::string( 1, str.sname ) ) != args_.end() ) {
        CG_INFO( "ArgumentsParser" ) << help_message();
        exit( 0 );
      }
    }
    //--- first loop over all required parameters
    for ( auto& par : required_params_ ) {
      const auto key = find( args_.begin(), args_.end(),
        par.name.empty() ? par.name : "--"+par.name );
      const auto skey = find( args_.begin(), args_.end(),
        "-"+std::string( 1, par.sname ) );
      if ( key == args_.end() && skey == args_.end() ) {
        throw CG_FATAL( "ArgumentsParser" )
          << help_message() << "\n\t"
          << "The following parameter was not set: '" << par.name << "'.";
      }
      const auto value = ( key != args_.end() ) ? key+1 : skey+1;
      if ( value == args_.end() )
        throw CG_FATAL( "ArgumentsParser" )
          << "Invalid value for parameter: " << par.name << ".";

      par.value = *value;
      par.parse();
    }
    //--- then loop over the optional parameters and assign the value
    //    for the variable with default value if not found from the
    //    arguments list
    for ( auto& par : optional_params_ ) {
      const auto key = find( args_.begin(), args_.end(),
        "--"+par.name );
      const auto skey = find( args_.begin(), args_.end(),
        "-"+std::string( 1, par.sname ) );
      if ( key != args_.end() || skey != args_.end() ) { // Parameter set
        const auto value = ( key != args_.end() ) ? key + 1 : skey + 1;
        if ( value != args_.end() ) {
          for ( const auto& par2 : optional_params_ )
            if ( *value == "--"+par2.name || *value == "-"+std::string( 1, par.sname ) )
              throw CG_FATAL( "ArgumentsParser" )
                << "Invalid value for parameter: " << par.name << ".";
          par.value = *value;
        }
        else if ( par.bool_variable )
          par.value = "1"; // if the flag is set, enabled by default
        else
          throw CG_FATAL( "ArgumentsParser" )
            << "Invalid value for parameter: " << par.name << ".";
      }
      par.parse();
    }
  }

  std::string
  ArgumentsParser::operator[]( std::string name ) const
  {
    for ( const auto& par : required_params_ ) {
      if ( "--"+par.name == name ) return par.value;
      if ( par.sname != '\0' && "-"+std::string( 1, par.sname ) == name )
        return par.value;
    }
    for ( const auto& par : optional_params_ ) {
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
    oss << "Usage: " << command_name_;
    for ( const auto& par : required_params_ ) {
      oss << ( !par.name.empty() ? " --" : " <arg>" ) << par.name;
      if ( par.sname != '\0' ) oss << "|-" << par.sname;
    }
    for ( const auto& par : optional_params_ ) {
      oss << " [--" << par.name;
      if ( par.sname != '\0' ) oss << " | -" << par.sname;
      oss << "]";
    }
    if ( required_params_.size() > 0 ) {
      oss << "\n    required argument" << ( ( required_params_.size() > 1 ) ? "s" : "" ) << ":";
      for ( const auto& par : required_params_ )
        oss << Form( ( par.sname != '\0' )
          ? "\n\t%s/-%1s\t%-28s"
          : "\n\t%s/%2s\t%-28s",
          ( !par.name.empty() ? "--"+par.name : "<arg>" ).c_str(),
          &par.sname, par.description.c_str() );
    }
    if ( optional_params_.size() > 0 ) {
      oss << "\n    optional argument" << ( ( optional_params_.size() > 1 ) ? "s" : "" ) << ":";
      for ( const auto& par : optional_params_ )
        oss << Form( ( par.sname != '\0' )
          ? "\n\t--%s/-%1s\t%-28s\tdefault = '%s'"
          : "\n\t--%s/%2s\t%-28s\tdefault = '%s'",
          par.name.c_str(), &par.sname, par.description.c_str(), par.value.c_str() );
    }
    oss << std::endl;
    return oss.str();
  }

  //----- simple parameters

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::string default_value, std::string* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( default_value ),
    str_variable( var ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::string* var, char sname ) :
    Parameter( name, description, "", var, sname )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, unsigned int default_value, unsigned int* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( std::to_string( default_value ) ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( var ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, unsigned int* var, char sname ) :
    Parameter( name, description, 0, var, sname )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, int default_value, int* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( Form( "%+i", default_value ) ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( var ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, int* var, char sname ) :
    Parameter( name, description, 0, var, sname )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, bool default_value, bool* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( Form( "%d", default_value ) ),
    str_variable( nullptr ), float_variable( nullptr ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( var ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, bool* var, char sname ) :
    Parameter( name, description, 0, var, sname )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, double default_value, double* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( Form( "%g", default_value ) ),
    str_variable( nullptr ), float_variable( var ),
    int_variable( nullptr ), uint_variable( nullptr ), bool_variable( nullptr ),
    vec_str_variable( nullptr ), vec_int_variable( nullptr ), vec_float_variable( nullptr )
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, double* var, char sname ) :
    Parameter( name, description, 0., var, sname )
  {}

  //----- vector of parameters

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<std::string> default_value, std::vector<std::string>* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( "" ),
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
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<int> default_value, std::vector<int>* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( "" ),
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
  {}

  ArgumentsParser::Parameter::Parameter( std::string name, std::string description, std::vector<double> default_value, std::vector<double>* var, char sname ) :
    name( name ), sname( sname ), description( description ),
    value( "" ),
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
  {}

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
                        && strcasecmp( "false", value.c_str() ) != 0 );
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

