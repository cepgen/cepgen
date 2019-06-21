#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"

using namespace std;

int main( int argc, char* argv[] )
{
  if ( argc < 2 )
    throw CG_FATAL( "main" ) << "Usage: " << argv[0] << " input-card";

  cepgen::Generator gen;
  gen.setParameters( cepgen::card::Handler::parse( argv[1] ) );

  cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debugInsideLoop;

  vector<double> x( 12, 0.3 );

  cout << "point: ";
  string delim;
  for ( const auto& v : x )
    cout << delim << v, delim = ", ";
  cout << endl;
  cout << "weight: " << gen.computePoint( &x[0] ) << endl;

  return 0;
}
