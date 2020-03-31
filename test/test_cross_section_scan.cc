#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Processes/Process.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/ArgumentsParser.h"

#include <fstream>

using namespace std;

int main( int argc, char* argv[] )
{
  string proc_name, output_file, scan;
  int proc_mode, npoints;
  double min_value, max_value, sqrts;
  double min_eta, max_eta;

  cepgen::ArgumentsParser parser( argc, argv );
  parser
    .addArgument( "proc-name", "name of the process to scan", &proc_name, 'p' )
    .addArgument( "min-value", "minimum value of scan", &min_value, 'l' )
    .addArgument( "max-value", "maximum value of scan", &max_value, 'H' )
    .addOptionalArgument( "proc-mode", "kin. mode to consider", 0, &proc_mode, 'm' )
    .addOptionalArgument( "scan", "type of scan to perform", "minpt", &scan, 's' )
    .addOptionalArgument( "num-points", "number of points to consider", 10, &npoints, 'n' )
    .addOptionalArgument( "output", "output file", "xsect.dat", &output_file, 'o' )
    .addOptionalArgument( "sqrts", "c.o.m. energy", 13.e3, &sqrts, 'w' )
    .addOptionalArgument( "min-eta", "minimum central system pseudorapidity", -2.5, &min_eta )
    .addOptionalArgument( "max-eta", "maximum central system pseudorapidity", +2.5, &max_eta )
    .parse();

  cepgen::Generator mg;

  //cepgen::Logger::get().level = cepgen::Logger::Level::error;

  cepgen::Parameters& par = mg.parameters();
  par.kinematics.setSqrtS( sqrts );
  par.kinematics.cuts.central.eta_single = { min_eta, max_eta };
  par.kinematics.cuts.remnants.mass_single.max() = 1000.0;
  par.setProcess( cepgen::proc::ProcessesFactory::get().build( proc_name ) );
  if ( proc_mode > 0 )
    par.kinematics.mode = static_cast<cepgen::KinematicsMode>( proc_mode );
  if ( !parser.extra_config().empty() )
    cepgen::card::CardsHandlerFactory::get().build( "cmd",
      cepgen::ParametersList().set<std::vector<std::string> >( "args", parser.extra_config() ) )->parse( "", par );
  CG_INFO( "main" ) << &par;

  double xsect, err_xsect;

  ofstream xsect_file( output_file );
  if ( !xsect_file.is_open() )
    throw CG_FATAL( "main" ) << "Output file \"" << output_file << "\" cannot be opened!";
  xsect_file << "# " << scan << "\txsect (pb)\td(xsect) (pb)\n";

  for ( int i = 0; i < npoints; i++ ) {
    const double value = min_value+( max_value-min_value )*i/npoints;
    if ( scan == "minpt" )
      par.kinematics.cuts.central.pt_single.min() = value;
    else if ( scan == "maxpt" )
      par.kinematics.cuts.central.pt_single.max() = value;
    CG_LOG( "main" ) << "Scan of \"" << scan << "\". Value = " << value << ".";
    mg.computeXsection( xsect, err_xsect );
    string out_line = cepgen::utils::format( "%.2f\t%.5f\t%.5f\n", value, xsect, err_xsect );
    xsect_file << out_line;
    cout << out_line;
    xsect_file.flush();
  }

  return 0;
}
