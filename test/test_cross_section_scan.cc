#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Cards/Handler.h"
#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Processes/Process.h"

#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/ArgumentsParser.h"

#include <fstream>

using namespace std;

int main( int argc, char* argv[] )
{
  string input_config, output_file, scan;
  int proc_mode, npoints;
  double min_value, max_value;

  cepgen::ArgumentsParser parser( argc, argv );
  parser
    .addArgument( "config", "base configuration", &input_config, 'i' )
    .addArgument( "min-value", "minimum value of scan", &min_value, 'l' )
    .addArgument( "max-value", "maximum value of scan", &max_value, 'H' )
    .addOptionalArgument( "proc-mode", "kin. mode to consider", 0, &proc_mode, 'm' )
    .addOptionalArgument( "scan", "type of scan to perform", "ptmin", &scan, 's' )
    .addOptionalArgument( "num-points", "number of points to consider", 10, &npoints, 'n' )
    .addOptionalArgument( "output", "output file", "xsect.dat", &output_file, 'o' )
    .parse();

  cepgen::Generator mg;
  mg.setParameters( cepgen::card::Handler::parse( input_config ) );

  cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::error;

  auto& par = mg.parameters();
  if ( proc_mode > 0 )
    par.kinematics.mode = static_cast<cepgen::KinematicsMode>( proc_mode );
  if ( !parser.extra_config().empty() )
    par = cepgen::card::CardsHandlerFactory::get().build( "cmd",
      cepgen::ParametersList().set<std::vector<std::string> >( "args", parser.extra_config() ) )->parse( "", par );
  CG_INFO( "main" ) << &par;

  double xsect, err_xsect;

  ofstream xsect_file( output_file );
  if ( !xsect_file.is_open() )
    throw CG_FATAL( "main" ) << "Output file \"" << output_file << "\" cannot be opened!";
  xsect_file << "# " << scan << "\txsect (pb)\td(xsect) (pb)\n";

  for ( int i = 0; i < npoints; i++ ) {
    const double value = min_value+( max_value-min_value )*i/npoints;
    if ( scan == "ptmin" )
      par.kinematics.cuts.central.pt_single.min() = value;
    else if ( scan == "ptmax" )
      par.kinematics.cuts.central.pt_single.max() = value;
    else if ( scan == "q2min" )
      par.kinematics.cuts.initial.q2.min() = value;
    else if ( scan == "q2max" )
      par.kinematics.cuts.initial.q2.max() = value;
    else if ( scan == "wmin" )
      par.kinematics.cuts.central.mass_sum.min() = value;
    else if ( scan == "wmax" )
      par.kinematics.cuts.central.mass_sum.max() = value;
    else if ( scan == "mxmin" )
      par.kinematics.cuts.remnants.mass_single.min() = value;
    else if ( scan == "mxmax" )
      par.kinematics.cuts.remnants.mass_single.max() = value;
    else if ( scan == "abseta" ) {
      par.kinematics.cuts.central.eta_single.min() = -value;
      par.kinematics.cuts.central.eta_single.max() = +value;
    }
    else if ( scan == "absrap" ) {
      par.kinematics.cuts.central.rapidity_single.min() = -value;
      par.kinematics.cuts.central.rapidity_single.max() = +value;
    }
    else if ( scan == "mpart" ) {
      auto prop = cepgen::PDG::get()( par.process().event()[cepgen::Particle::CentralSystem][0].pdgId() );
      prop.mass = value;
      cepgen::PDG::get().define( prop );
      par.process().clear();
    }
    else
      throw CG_FATAL( "main" ) << "Invalid variable to be scanned: \"" << scan << "\"!";
    CG_LOG( "main" ) << "Scan of \"" << scan << "\". Value = " << value << ".";
    mg.computeXsection( xsect, err_xsect );
    string out_line = cepgen::utils::format( "%.2f\t%.8e\t%.8e\n", value, xsect, err_xsect );
    xsect_file << out_line;
    cout << out_line;
    xsect_file.flush();
  }

  return 0;
}
