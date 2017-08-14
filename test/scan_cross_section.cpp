#include "CepGen/Generator.h"

using namespace std;

int main( int argc, char* argv[] )
{
  if ( argc<5 ) {
    InError( Form( "Usage: %s <process mode=1..4> <num points> <min value> <max value> [output file=xsect.dat]", argv[0] ) );
    return -1;
  }
  const unsigned int proc_mode = atoi( argv[1] ),
                     npoints = atoi( argv[2] );
  const float min_value = atof( argv[3] ),
              max_value = atof( argv[4] );
  const char* output_file = ( argc>5 ) ? argv[5] : "xsect.dat";

  CepGen::Generator mg;

  CepGen::Logger::get().level = CepGen::Logger::Error;

  CepGen::Parameters* par = mg.parameters.get();
  par->kinematics.eta_min = -2.5; par->kinematics.eta_max = 2.5;
  par->kinematics.in1p = par->kinematics.in2p = 6.5e3;
  par->kinematics.mx_max = 1000.0;
  par->setProcess( new CepGen::Process::GamGamLL );
  par->process_mode = static_cast<CepGen::Kinematics::ProcessMode>( proc_mode );
  par->dump();

  double xsect, err_xsect;

  ofstream xsect_file( output_file );
  if ( !xsect_file.is_open() ) {
    InError( Form( "Output file \"%s\" cannot be opened!", output_file ) );
    return -2;
  }

  for ( unsigned int i=0; i<npoints; i++ ) {
    par->kinematics.pt_min = min_value + (max_value-min_value)*i/npoints;
    mg.computeXsection( xsect, err_xsect );
    xsect_file << Form( "%.2f\t%.5f\t%.5f\n", par->kinematics.pt_min, xsect, err_xsect );
    xsect_file.flush();
  }

  return 0;
}
