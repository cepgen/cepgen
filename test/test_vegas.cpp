#include "CepGen/Generator.h"
#include "CepGen/Processes/TestProcess.h"

#include <assert.h>

int
main( int argc, char* argv[] )
{
  const double num_sigma = 3.0;
  const double exact = 1.3932039296856768591842462603255;

  Timer tmr;
  CepGen::Generator mg;

  //CepGen::Logger::get().level = CepGen::Logger::Debug;

  mg.parameters->setProcess( new CepGen::Process::TestProcess );
  //mg.parameters->vegas.ncvg = 50000;
  //mg.parameters->vegas.itvg = 5;

  Information( Form( "Initial configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
  tmr.reset();

  double result, error;
  mg.computeXsection( result, error );

  //const double sigma = ( fabs( xsec_ref-xsec_cepgen ) ) / sqrt( err_xsec_cepgen*err_xsec_cepgen + err_xsec_ref*err_xsec_ref );

  //Information( Form( "Computed cross section:\n\tRef.   = %.3e +/- %.3e\n\tCepGen = %.3e +/- %.3e\n\tPull: %.6f", xsec_ref, err_xsec_ref, xsec_cepgen, err_xsec_cepgen, sigma ) );

  //Information( Form( "Computation time: %.3f ms", tmr.elapsed()*1.e3 ) );
  tmr.reset();

  assert( fabs( exact - result ) < num_sigma * error );

  Information( "ALL TESTS PASSED!" );

  return 0;
}
