#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main( int argc, char* argv[] )
{
  string input_card, point;
  bool debug;

  cepgen::ArgumentsParser( argc, argv )
    .addArgument( "input", "input card", &input_card, 'i' )
    .addOptionalArgument( "point", "point to test", "", &point, 'p' )
    .addOptionalArgument( "debug", "debugging mode", false, &debug, 'd' )
    .parse();

  if ( !debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debugInsideLoop;

  vector<double> x;
  if ( !point.empty() ) {
    stringstream iss( point );
    double buf;
    while ( iss >> buf )
      x.emplace_back( buf );
    if ( x.size() < 2 )
      x = vector<double>( 12, x[0] );
    x.resize( 12 );
  }
  else
    x = vector<double>( 12, 0.3 );

  cepgen::Generator gen;
  gen.setParameters( cepgen::card::Handler::parse( input_card.c_str() ) );
  CG_INFO( "main" ) << gen.parametersPtr();

  cout << "point: ";
  string delim;
  for ( const auto& v : x )
    cout << delim << v, delim = ", ";
  cout << endl;
  cout << "weight: " << gen.computePoint( &x[0] ) << endl;

  return 0;
}
