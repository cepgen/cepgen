#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Processes/GamGamLL.h"

#include <fstream>

using namespace std;

int main( int argc, char* argv[] )
{
  if ( argc < 5 )
    throw CG_FATAL( "main" )
      << "Usage:\n"
      << argv[0] << " <process mode=1..4> <num points> <min value> <max value> [output file=xsect.dat]";

  const unsigned int proc_mode = atoi( argv[1] ),
                     npoints = atoi( argv[2] );
  const float min_value = atof( argv[3] ),
              max_value = atof( argv[4] );
  const char* output_file = ( argc>5 ) ? argv[5] : "xsect.dat";

  CepGen::Generator mg;

  //CepGen::Logger::get().level = CepGen::Logger::Level::error;

  CepGen::Parameters* par = mg.parameters.get();
  par->kinematics.cuts.central.eta_single = { -2.5, 2.5 };
  par->kinematics.cuts.remnants.mass_single.max() = 1000.0;
  par->setProcess( new CepGen::Process::GamGamLL );
  par->kinematics.mode = static_cast<CepGen::Kinematics::Mode>( proc_mode );
  par->dump();

  double xsect, err_xsect;

  ofstream xsect_file( output_file );
  if ( !xsect_file.is_open() )
    throw CG_FATAL( "main" ) << "Output file \"" << output_file << "\" cannot be opened!";

  for ( unsigned int i=0; i<npoints; i++ ) {
    par->kinematics.cuts.central.pt_single.min() = min_value + (max_value-min_value)*i/npoints;
    //par->kinematics.dump();
    mg.computeXsection( xsect, err_xsect );
    xsect_file << CepGen::Form( "%.2f\t%.5f\t%.5f\n", par->kinematics.cuts.central.pt_single.min(), xsect, err_xsect );
    cout << CepGen::Form( "%.2f\t%.5f\t%.5f\n", par->kinematics.cuts.central.pt_single.min(), xsect, err_xsect );
    xsect_file.flush();
  }

  return 0;
}
