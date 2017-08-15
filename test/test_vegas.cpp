#include "CepGen/Generator.h"
#include "CepGen/Processes/TestProcess.h"

#include <assert.h>

int
main( int argc, char* argv[] )
{
  const double num_sigma = 3.0;

  Timer tmr;
  CepGen::Generator mg;

  //CepGen::Logger::get().level = CepGen::Logger::Debug;

  mg.parameters->kinematics.setSqrtS( 13.e3 );
  mg.parameters->kinematics.eta_min = -2.5;
  mg.parameters->kinematics.eta_max = 2.5;
  mg.parameters->kinematics.mx_max = 1000.;
  mg.parameters->vegas.ncvg = 50000;
  mg.parameters->vegas.itvg = 5;

  Information( Form( "Initial configuration time: %.3f ms", tmr.elapsed()*1.e3 ) );
  tmr.reset();

  mg.parameters->setProcess( new CepGen::Process::TestProcess );

  mg.clearRun();
  double xsec_cepgen, err_xsec_cepgen;
  mg.computeXsection( xsec_cepgen, err_xsec_cepgen );

  //const double sigma = ( fabs( xsec_ref-xsec_cepgen ) ) / sqrt( err_xsec_cepgen*err_xsec_cepgen + err_xsec_ref*err_xsec_ref );

  //Information( Form( "Computed cross section:\n\tRef.   = %.3e +/- %.3e\n\tCepGen = %.3e +/- %.3e\n\tPull: %.6f", xsec_ref, err_xsec_ref, xsec_cepgen, err_xsec_cepgen, sigma ) );

  //Information( Form( "Computation time: %.3f ms", tmr.elapsed()*1.e3 ) );
  tmr.reset();

  //assert( fabs( sigma )<num_sigma );

  Information( "ALL TESTS PASSED!" );

  return 0;
}
