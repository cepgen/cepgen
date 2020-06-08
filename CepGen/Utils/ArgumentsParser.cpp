#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Core/Exception.h"

#include <sstream>
#include <algorithm>
#include <cstring>

namespace cepgen
{
  ArgumentsParser::ArgumentsParser( int argc, char* argv[] ) :
    help_str_( { { "help", 'h' } } ),
    config_str_( { { "cmd", 'c' } } )
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
    std::vector<std::pair<std::string,std::string> > arguments;
    if ( !args_.empty() ) {
      //--- check if help message is requested
      for ( const auto& str : help_str_ )
        if ( find( args_.begin(), args_.end(), "--"+str.name ) != args_.end()
          || find( args_.begin(), args_.end(), "-"+std::string( 1, str.sname ) ) != args_.end() ) {
          print_help();
          exit(0);
        }
      for ( const auto& str : config_str_ ) {
        auto cfg_long = find( args_.begin(), args_.end(), "--"+str.name );
        if ( cfg_long != args_.end() )
          extra_config_ = std::vector<std::string>( cfg_long+1, args_.end() );
        auto cfg_short = find( args_.begin(), args_.end(), "-"+std::string( 1, str.sname ) );
        if ( cfg_short != args_.end() )
          extra_config_ = std::vector<std::string>( cfg_short+1, args_.end() );
      }
      if ( !extra_config_.empty() )
        args_.resize( args_.size()-extra_config_.size()-1 );
      for ( auto it_arg = args_.begin(); it_arg != args_.end(); ++it_arg ) {
        auto arg_val = utils::split( *it_arg, '=' ); // particular case for --arg=value
        if ( arg_val.size() == 1 )
          arg_val.emplace_back( it_arg != std::prev( args_.end() ) : *std::next( it_arg ) : "" );
        arguments.emplace_back( std::make_pair( arg_val.at( 0 ), arg_val.at( 1 ) ) );
      }
      CG_WARNING("")<<arguments;
    }
    //--- loop over all parameters
    size_t i = 0;
    for ( auto& par : params_ ) {
      if ( par.name.empty() ) {
        //--- no argument name ; fetching by index
        if ( i >= args_.size() )
          throw CG_FATAL( "ArgumentsParser" )
            << help_message()
            << " Failed to retrieve required <arg" << i << ">.";
        par.value = !par.bool_variable ? args_.at( i ) : "1";
      }
      else {
        auto it_key = find( args_.begin(), args_.end(), "--"+par.name );
        if ( it_key == args_.end() ) // allow for '-param' instead of '--param'
          it_key = find( args_.begin(), args_.end(), "-"+par.name );
        const auto it_skey = find( args_.begin(), args_.end(), "-"+std::string( 1, par.sname ) );
        CG_WARNING("")<<par.name<<":"<<args_<<":"<<(it_key!=args_.end())<<"|"<<(it_skey!=args_.end());
        if ( it_key != args_.end() || it_skey != args_.end() ) {
          // matched a long/short argument
          const auto it_value = it_key != args_.end()
            ? std::next( it_key )
            : std::next( it_skey );
          CG_WARNING("")<<par<<":"<<*it_key<<"="<<*it_value;
          if ( it_value != args_.end() )
            par.value = *it_value;
          else if ( par.bool_variable )
            par.value = "1"; // if the flag is set, enabled by default
          else
            throw CG_FATAL( "ArgumentsParser" )
              << "Invalid value for parameter: '" << par.name << "'.";
          ++i;
        }
        else if ( args_.size() > i && args_.at( i )[0] != '-' )
          par.value = args_.at( i );
        else if ( !par.optional ) // no match
          throw CG_FATAL( "ArgumentsParser" )
            << help_message()
            << " The following parameter was not set: '" << par.name << "'.";
      }
      CG_INFO( "ArgumentsParser" )
        << "Parameter '" << i << "|--" << par.name << "|-" << par.sname << "'"
        << " has value '" << par.value << "'.";
      par.parse();
      ++i;
    }
    /*auto it_param = params_.begin();
    for ( auto it_arg = args_.begin(); it_arg != args_.end(); ++it_arg ) {
      auto arg_val = utils::split( *it_arg, '=' ); // particular case for --arg=value
      if ( arg_val.size() == 1 && it_arg != std::prev( args_.end() ) )
        arg_val.emplace_back( *std::next( it_arg ) );
      auto param = parameter( arg_val.at( 0 ) );
      if ( !param )
        param = &(*it_param++);
      if ( !param )
        throw CG_FATAL( "ArgumentsParser" )
          << "No correspondance found for argument \"" << *it_arg << "\"!";
      param->value = arg_val.at( 1 );
      if ( param->bool_variable ) { // particular case for boolean flags
        const auto word = utils::tolower( arg_val.at( 1 ) );
        if ( word == "on" || word != "off"
          || word == "yes" || word != "no"
          || word == "true" || word != "false" ) {
          param->value = "1";
          ++it_arg;
        }
      }
      else
        ++it_arg;
      param->parse();
      CG_WARNING("")<<arg_val<<"-->"<<*param;
    }*/
    return *this;
  }

  ArgumentsParser::Parameter*
  ArgumentsParser::parameter( const std::string& arg )
  {
    for ( auto& par : params_ ) {
      if ( arg.find( "--" ) != arg.npos && arg == "--"+par.name )
        return &par;
      if ( arg.find( "-" ) != arg.npos && arg == "-"+std::string( 1, par.sname ) )
        return &par;
    }
    return nullptr;
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
        oss << " [";
      }
      else {
        req_params.emplace_back( std::make_pair( par, i ) );
        oss << " ";
      }
      oss << ( !par.name.empty() ? "--" : " <arg"+std::to_string( i )+">" ) << par.name;
      if ( par.sname != '\0' )
        oss << ( !par.name.empty() ? "|" : "" ) << "-" << par.sname;
      if ( par.optional )
        oss << "]";
      ++i;
    }
    if ( req_params.size() > 0 ) {
      oss << "\n    " << utils::s( "required argument", req_params.size(), false ) << ":";
      for ( const auto& par : req_params )
        oss << utils::format( "\n\t%s%-18s\t%-28s",
          ( par.first.sname != '\0' ? "-"+std::string( 1, par.first.sname )+"|" : "" ).c_str(),
          ( !par.first.name.empty() ? "--"+par.first.name : "<arg"+std::to_string( par.second )+">" ).c_str(),
          par.first.description.c_str() );
    }
    if ( opt_params.size() > 0 ) {
      oss << "\n    " << utils::s( "optional argument", opt_params.size(), false ) << ":";
      for ( const auto& par : opt_params )
        oss << utils::format( "\n\t%s%-18s\t%-28s\tdefault = '%s'",
          ( par.first.sname != '\0' ? "-"+std::string( 1, par.first.sname )+"|" : "" ).c_str(),
          ( !par.first.name.empty() ? "--"+par.first.name : "<arg"+std::to_string( par.second )+">" ).c_str(),
          par.first.description.c_str(), par.first.value.c_str() );
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
    value( utils::format( "%+i", default_value ) ), optional( true ),
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
    value( utils::format( "%d", default_value ) ), optional( true ),
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
    value( utils::format( "%g", default_value ) ), optional( true ),
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
      value += ( ( ( i++ > 0 ) ? "," : "" )+utils::format( "%d", var ) );
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
      value += ( ( ( i++ > 0 ) ? "," : "" )+utils::format( "%g", flt ) );
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
      try {
        *float_variable = std::stod( value );
      } catch ( const std::invalid_argument& ) {
        throw CG_FATAL( "ArgumentsParser:Parameter:parse" )
          << "Failed to parse variable '" << name << "' as float!";
      }
    else if ( int_variable != nullptr )
      try {
        *int_variable = std::stoi( value );
      } catch ( const std::invalid_argument& ) {
        throw CG_FATAL( "ArgumentsParser:Parameter:parse" )
          << "Failed to parse variable '" << name << "' as integer!";
      }
    else if ( uint_variable != nullptr )
      try {
        *uint_variable = std::stoi( value );
      } catch ( const std::invalid_argument& ) {
        throw CG_FATAL( "ArgumentsParser:Parameter:parse" )
          << "Failed to parse variable '" << name << "' as unsigned integer!";
      }
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

  std::ostream&
  operator<<( std::ostream& os, const ArgumentsParser::Parameter& par )
  {
    return os << "Parameter{"
      << "--" << par.name
      << ",-" << std::string( 1, par.sname )
      << ( !par.description.empty() ? ","+par.description : "" )
      << par.value
      << ",opt:" << std::boolalpha << par.optional
      << "}";
  }
}

