#include "CepGen/Core/Functional.h"
#include "CepGen/Core/Exception.h"

#include <math.h>
#include <string>
#include <iostream>

using namespace std;

int
main()
{
  const double epsilon = 1.e-9; // tolerance

  { // test with a 1-variable function
    const double exp_result_test1 = 6.795704571;
    cepgen::utils::Functional<1> test( "2.5*exp(0.1*x)", { "x" } );
    if ( fabs( test.eval( 10. ) - exp_result_test1 ) > epsilon  )
      throw CG_FATAL( "main" ) << "Test 1.1 failed!";
    if ( fabs( test.eval( { 10. } ) - exp_result_test1 ) > epsilon  )
      throw CG_FATAL( "main" ) << "Test 1.2 failed!";
    cout << "Test 1 passed!" << endl;
  }
  { // test with an invalid function
    bool passed = false;
    try {
      cepgen::utils::Functional<1> test( "sqrt(x+x**3-log(10)", { "x" } );
      test.eval( 10 );
    } catch ( const cepgen::Exception& e ) { cout << "Test 2 passed!" /*<< e.what()*/ << endl; passed = true; }
    if ( !passed )
      throw CG_FATAL( "main" ) << "Test 2 failed!";
  }
  { // test with a 2-variables function
    cepgen::utils::Functional<2> test( "sqrt(a^2+b^2)", { "a", "b" } );
    if ( fabs( test.eval( { 3, 4 } ) - 5.0 ) > epsilon  )
      throw CG_FATAL( "main" ) << "Test 3 failed!";
    cout << "Test 3 passed!" << endl;
  }
  { // test with an invalid function
    cepgen::utils::Functional<1> test( "a**2", { "a" } );
    bool passed = true;
    try { test.eval( 10 ); passed = false; } catch ( const cepgen::Exception& e ) {}
    try { test.eval( { 10 } ); passed = false; } catch ( const cepgen::Exception& e ) {}
    if ( !passed  )
      throw CG_FATAL( "main" ) << "Test 4 failed!";
    cout << "Test 4 passed!" << endl;
  }

  return 0;
}
