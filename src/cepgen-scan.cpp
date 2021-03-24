#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"

#include "CepGen/Core/ExportModule.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Cards/Handler.h"
#include "CepGen/Processes/Process.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Modules/CardsHandlerFactory.h"

#include "CepGen/Utils/Logger.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/ArgumentsParser.h"

#include <fstream>

using namespace std;

int main( int argc, char* argv[] )
{
  string input_config, output_file, scan;
  int npoints;
  double min_value, max_value;
  vector<double> points;
  bool debug;

  cepgen::ArgumentsParser parser( argc, argv );
  parser
    .addArgument( "config,i", "base configuration", &input_config )
    .addOptionalArgument( "scan,s", "type of scan to perform", &scan, "ptmin" )
    .addOptionalArgument( "min,l", "minimum value of scan", &min_value, 1. )
    .addOptionalArgument( "max,H", "maximum value of scan", &max_value, 11. )
    .addOptionalArgument( "num-points,n", "number of points to consider", &npoints, 10 )
    .addOptionalArgument( "points,p", "list of points to consider", &points, vector<double>{} )
    .addOptionalArgument( "output,o", "output file", &output_file, "xsect.dat" )
    .addOptionalArgument( "debug,d", "debugging mode", &debug, false )
    .parse();

  if ( debug )
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debug;
  else
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::nothing;

  cepgen::Generator mg;
  mg.setParameters( cepgen::card::Handler::parse( input_config ) );

  if ( !parser.extra_config().empty() )
    mg.setParameters( cepgen::card::CardsHandlerFactory::get().build( "cmd",
      cepgen::ParametersList().set<std::vector<std::string> >( "args", parser.extra_config() ) )
      ->parse( "", mg.parametersPtr() ) );

  CG_INFO( "main" ) << mg.parameters();

  double xsect, err_xsect;

  ofstream xsect_file( output_file );
  if ( !xsect_file.is_open() )
    throw CG_FATAL( "main" ) << "Output file \"" << output_file << "\" cannot be opened!";
  xsect_file << "# " << scan << "\txsect (pb)\td(xsect) (pb)\n";

  auto& par = mg.parametersRef();
  //--- ensure nothing is written in the output sequence
  par.outputModulesSequence().clear();

  if ( points.empty() )
    for ( int i = 0; i <= npoints; ++i )
      points.emplace_back( min_value+( max_value-min_value )*i/npoints );

  for ( const auto& value : points ) {
    if ( scan == "ptmin" )
      par.kinematics.cuts.central.pt_single().min() = value;
    else if ( scan == "ptmax" )
      par.kinematics.cuts.central.pt_single().max() = value;
    else if ( scan == "q2min" )
      par.kinematics.cuts.initial.q2().min() = value;
    else if ( scan == "q2max" )
      par.kinematics.cuts.initial.q2().max() = value;
    else if ( scan == "wmin" )
      par.kinematics.cuts.central.mass_sum().min() = value;
    else if ( scan == "wmax" )
      par.kinematics.cuts.central.mass_sum().max() = value;
    else if ( scan == "mxmin" )
      par.kinematics.cuts.remnants.mx().min() = value;
    else if ( scan == "mxmax" )
      par.kinematics.cuts.remnants.mx().max() = value;
    else if ( scan == "abseta" ) {
      par.kinematics.cuts.central.eta_single().min() = -value;
      par.kinematics.cuts.central.eta_single().max() = +value;
    }
    else if ( scan == "absrap" ) {
      par.kinematics.cuts.central.rapidity_single().min() = -value;
      par.kinematics.cuts.central.rapidity_single().max() = +value;
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
