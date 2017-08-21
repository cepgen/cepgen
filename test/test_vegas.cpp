#include "CepGen/Generator.h"
#include "CepGen/Processes/TestProcess.h"

#include <iostream>

using namespace std;

int
main( int argc, char* argv[] )
{
  const double exact = 1.3932039296856768591842462603255;

  CepGen::Generator mg;

  //CepGen::Logger::get().level = CepGen::Logger::Debug;

  mg.parameters->setProcess( new CepGen::Process::TestProcess );
  //mg.parameters->vegas.ncvg = 50000;
  //mg.parameters->vegas.itvg = 5;

  double result, error;
  mg.computeXsection( result, error );

  if ( fabs( exact - result ) > 1.0 * error ) return -1;

  cout << "Test 1 passed!" << endl;

  return 0;
}
