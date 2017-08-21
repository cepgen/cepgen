#include "CepGen/Cards/FunctionBuilder.h"

#include <string>
#include <iostream>

using namespace std;

int
main()
{
  const double epsilon = 1.0e-6;
  {
    const double exp_result_test1 = 6.795704571;
    CepGen::FunctionBuilder<1> test1( "2.5*exp(0.1*x)", { "x" } );
    if ( fabs( test1.eval( 10. ) - exp_result_test1 ) > epsilon ) return -1;
    cout << "Test 1 passed!" << endl;
  }
  {
    const double exp_result_test2 = 5.0;
    CepGen::FunctionBuilder<2> test2( "sqrt(a^2+b^2)", { "a", "b" } );
    if ( fabs( test2.eval( { 3, 4 } ) - exp_result_test2 ) > epsilon ) return -1;
    cout << "Test 2 passed!" << endl;
  }

  return 0;
}
