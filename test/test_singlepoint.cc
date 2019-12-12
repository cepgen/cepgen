#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main( int argc, char* argv[] )
{
  const size_t ps_size = 12;
  string input_card;
  vector<double> point;
  bool debug;

  cepgen::ArgumentsParser( argc, argv )
    .addArgument( "input", "input card", &input_card, 'i' )
    .addOptionalArgument( "point", "point to test", vector<double>( ps_size, 0.3 ), &point, 'p' )
    .addOptionalArgument( "debug", "debugging mode", false, &debug, 'd' )
    .parse();

  if ( point.size() < 2 ) {
    point = vector<double>( ps_size, point[0] );
    point.resize( ps_size );
  }

  cepgen::Generator gen;
  gen.setParameters( cepgen::card::Handler::parse( input_card )->parameters() );
  CG_INFO( "main" ) << gen.parametersPtr();

  if ( !debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debugInsideLoop;

  cout << "point: ";
  string delim;
  for ( const auto& v : point )
    cout << delim << v, delim = ", ";
  cout << endl;
  cout << "weight: " << gen.computePoint( &point[0] ) << endl;

  return 0;
}
