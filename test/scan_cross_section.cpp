#include "core/MCGen.h"

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

  MCGen mg;

  Logger::GetInstance()->Level = Logger::Error;

  Parameters* par = mg.parameters;
  par->mineta = -2.5; par->maxeta = 2.5;
  par->in1p = par->in2p = 6.5e3;
  par->maxmx = 1000.0;
  par->process = new GamGamLL;
  par->process_mode = static_cast<Kinematics::ProcessMode>( proc_mode );
  par->Dump();

  double xsect, err_xsect;

  ofstream xsect_file( output_file );
  if ( !xsect_file.is_open() ) {
    InError( Form( "Output file \"%s\" cannot be opened!", output_file ) );
    return -2;
  }

  for ( unsigned int i=0; i<npoints; i++ ) {
    par->minpt = min_value + (max_value-min_value)*i/npoints;
    mg.ComputeXsection( &xsect, &err_xsect );
    xsect_file << Form( "%.2f\t%.5f\t%.5f\n", par->minpt, xsect, err_xsect );
    xsect_file.flush();
  }

  return 0;
}
