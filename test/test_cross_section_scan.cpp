#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Processes/ProcessesHandler.h"

#include <fstream>

using namespace std;

int main( int argc, char* argv[] )
{
  if ( argc < 6 )
    throw CG_FATAL( "main" )
      << "Usage:\n"
      << argv[0] << " <process name> <process mode=1..4> <num points> <min value> <max value> [output file=xsect.dat]";

  const string proc_name = argv[1];
  const unsigned int proc_mode = atoi( argv[2] ), npoints = atoi( argv[3] );
  const float min_value = atof( argv[4] ), max_value = atof( argv[5] );
  const char* output_file = ( argc > 6 ) ? argv[6] : "xsect.dat";

  cepgen::Generator mg;

  //cepgen::Logger::get().level = cepgen::Logger::Level::error;

  cepgen::Parameters& par = mg.parameters();
  par.kinematics.cuts.central.eta_single = { -2.5, 2.5 };
  par.kinematics.cuts.remnants.mass_single.max() = 1000.0;
  par.setProcess( cepgen::proc::ProcessesHandler::get().build( proc_name ) );
  par.kinematics.mode = static_cast<cepgen::KinematicsMode>( proc_mode );
  CG_INFO( "main" ) << &par;

  double xsect, err_xsect;

  ofstream xsect_file( output_file );
  if ( !xsect_file.is_open() )
    throw CG_FATAL( "main" ) << "Output file \"" << output_file << "\" cannot be opened!";

  for ( unsigned short i=0; i < npoints; i++ ) {
    par.kinematics.cuts.central.pt_single.min() = min_value+( max_value-min_value )*i/npoints;
    cout << par.kinematics.cuts.central.pt_single << endl;
    mg.computeXsection( xsect, err_xsect );
    xsect_file << cepgen::Form( "%.2f\t%.5f\t%.5f\n", par.kinematics.cuts.central.pt_single.min(), xsect, err_xsect );
    cout << cepgen::Form( "%.2f\t%.5f\t%.5f\n", par.kinematics.cuts.central.pt_single.min(), xsect, err_xsect );
    xsect_file.flush();
  }

  return 0;
}
