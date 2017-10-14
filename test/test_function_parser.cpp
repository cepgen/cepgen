#include "CepGen/Core/Functional.h"
#include "CepGen/Core/Exception.h"

#include <string>
#include <iostream>

using namespace std;

int
main()
{
  const double epsilon = 1.e-9; // tolerance

  { // test with a 1-variable function
    const double exp_result_test1 = 6.795704571;
    CepGen::Functional<1> test1( "2.5*exp(0.1*x)", { "x" } );
    if ( fabs( test1.eval( 10. ) - exp_result_test1 ) > epsilon  ) throw CepGen::Exception( __PRETTY_FUNCTION__, "Test 1.1 failed!", CepGen::FatalError );
    if ( fabs( test1.eval( { 10. } ) - exp_result_test1 ) > epsilon  ) throw CepGen::Exception( __PRETTY_FUNCTION__, "Test 1.2 failed!", CepGen::FatalError );
    cout << "Test 1 passed!" << endl;
  }
  { // test with a 2-variables function
    CepGen::Functional<2> test2( "sqrt(a^2+b^2)", { "a", "b" } );
    if ( fabs( test2.eval( { 3, 4 } ) - 5.0 ) > epsilon  ) throw CepGen::Exception( __PRETTY_FUNCTION__, "Test 2 failed!", CepGen::FatalError );
    cout << "Test 2 passed!" << endl;
  }
  { // test with an invalid function
    CepGen::Functional<1> test3( "a**2", { "a" } );
    bool passed = true;
    try { test3.eval( 10 ); passed = false; } catch ( CepGen::Exception& e ) {}
    try { test3.eval( { 10 } ); passed = false; } catch ( CepGen::Exception& e ) {}
    if ( !passed  ) throw CepGen::Exception( __PRETTY_FUNCTION__, "Test 3 failed!", CepGen::FatalError );
    cout << "Test 3 passed!" << endl;
  }

  return 0;
}
