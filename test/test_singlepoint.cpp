#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"

using namespace std;

int main( int argc, char* argv[] )
{
  if ( argc < 2 )
    throw CG_FATAL( "main" ) << "Usage: " << argv[0] << " input-card";

  vector<double> x;
  if ( argc > 2 ) {
    stringstream iss( argv[2] );
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
  gen.setParameters( cepgen::card::Handler::parse( argv[1] ) );
  CG_INFO( "main" ) << gen.parametersPtr();

  cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debugInsideLoop;

  cout << "point: ";
  string delim;
  for ( const auto& v : x )
    cout << delim << v, delim = ", ";
  cout << endl;
  cout << "weight: " << gen.computePoint( &x[0] ) << endl;

  return 0;
}
