#include "CepGen/Generator.h"
#include "CepGen/Processes/TestProcess.h"

#include <iostream>
#include <assert.h>

using namespace std;

int
main( int argc, char* argv[] )
{
  const double exact = 1.3932039296856768591842462603255;
  CepGen::Logger::get().level = CepGen::Logger::Nothing;

  CepGen::Generator mg;

  mg.parameters->setProcess( new CepGen::Process::TestProcess );
  mg.parameters->vegas.ncvg = 500000;
  //mg.parameters->vegas.itvg = 5;

  double result, error;
  mg.computeXsection( result, error );

  assert( fabs( exact - result ) < 2.0 * error );

  cout << "Test 1 passed!" << endl;

  return 0;
}
