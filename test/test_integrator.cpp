#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Core/Integrator.h"
#include "TestProcess.h"

#include <iostream>

using namespace std;

int
main( int argc, char* argv[] )
{
  if ( argc < 3 || string( argv[2] ) != "debug" )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;

  const double max_sigma = 3.0;

  cepgen::Generator mg;
  if ( argc > 1 && string( argv[1] ) == "plain" )
    mg.parameters->integrator.type = cepgen::IntegratorType::plain;
  if ( argc > 1 && string( argv[1] ) == "vegas" )
    mg.parameters->integrator.type = cepgen::IntegratorType::Vegas;
  if ( argc > 1 && string( argv[1] ) == "miser" )
    mg.parameters->integrator.type = cepgen::IntegratorType::MISER;

  double result, error;

  { // test 1
    const double exact = 1.3932039296856768591842462603255;
    mg.parameters->setProcess( new cepgen::proc::TestProcess<3> );
    mg.computeXsection( result, error );
    if ( fabs( exact - result ) > max_sigma * error )
      throw CG_FATAL( "main" ) << "pull = " << fabs( exact-result )/error << ".";
    cout << "Test 1 passed!" << endl;
  }
  { // test 2
    const double exact = 2./3.;
    mg.parameters->setProcess( new cepgen::proc::TestProcess<2>( "x^2+y^2", { "x", "y" } ) );
    mg.computeXsection( result, error );
    if ( fabs( exact - result ) > max_sigma * error )
      throw CG_FATAL( "main" ) << "pull = " << fabs( exact-result )/error << ".";
    cout << "Test 2 passed!" << endl;
  }
  { // test 3
    const double exact = 13./12.;
    mg.parameters->setProcess( new cepgen::proc::TestProcess<3>( "x+y^2+z^3", { "x", "y", "z" } ) );
    mg.computeXsection( result, error );
    if ( fabs( exact - result ) > max_sigma * error )
      throw CG_FATAL( "main" ) << "pull = " << fabs( exact-result )/error << ".";
    cout << "Test 3 passed!" << endl;
  }

  return 0;
}
