#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/FormFactors/Parameterisation.h"

#include "CepGen/Utils/ArgumentsParser.h"

#include <fstream>

#include "Canvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"

using namespace std;

int main( int argc, char* argv[] )
{
  int mode, strfun_type, num_points;
  double mx, q2min, q2max;
  string output_file;

  cepgen::ArgumentsParser( argc, argv )
    .addArgument( "mode,t", "beam modelling", &mode, (int)cepgen::mode::Beam::ProtonElastic )
    .addOptionalArgument( "mx,M", "diffractive mass (GeV/c^2)", &mx, 100. )
    .addOptionalArgument( "sf,s", "structure functions modelling", &strfun_type, 301 )
    .addOptionalArgument( "q2min,m", "minimal parton virtuality (GeV^2)", &q2min, 1. )
    .addOptionalArgument( "q2max,M", "maximal parton virtuality (GeV^2)", &q2max, 10000. )
    .addOptionalArgument( "npoints,n", "number of x-points to scan", &num_points, 500 )
    .addOptionalArgument( "output,o", "output file name", &output_file, "formfacs.scan.output.txt" )
    .parse();

  cepgen::initialise();
  //const double mi = cepgen::PDG::get().mass( cepgen::PDG::proton );
  //const double mi2 = mi*mi, mx2 = mx*mx;

  ofstream out( output_file );
  out
    << "# form factors: ";
  string sep;
  for ( const auto& fftype : cepgen::formfac::FormFactorsFactory::get().modules() )
    out << sep << fftype, sep = ", ";

  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build( strfun_type );
  out
    << "\n"
    << "# structure functions: " << sf.get() << "\n"
    << "# q2 in [" << q2min << ", " << q2max << "] GeV^2\n";

  vector<unique_ptr<cepgen::formfac::Parameterisation> > formfacs;
  vector<TGraph*> g_formfacs_fe, g_formfacs_fm;
  for ( const auto& fftype : cepgen::formfac::FormFactorsFactory::get().modules() ) {
    formfacs.emplace_back( cepgen::formfac::FormFactorsFactory::get().build( fftype ) );
    ( *formfacs.rbegin() )->setStructureFunctions( sf.get() );
    g_formfacs_fe.emplace_back( new TGraph );
    ( *g_formfacs_fe.rbegin() )->SetTitle( ( fftype+";Q^{2} (GeV^{2});F_{E}" ).c_str() );
    g_formfacs_fm.emplace_back( new TGraph );
    ( *g_formfacs_fm.rbegin() )->SetTitle( ( fftype+";Q^{2} (GeV^{2});F_{M}" ).c_str() );
  }
  for ( int i = 0; i < num_points; ++i ) {
    const double q2 = q2min+i*( q2max-q2min )/( num_points-1 );
    out << q2 << "\t";
    size_t j = 0;
    for ( auto& ff : formfacs ) {
      const auto formfac = ( *ff )( (cepgen::mode::Beam)mode, q2, mx );
      out << "\t" << formfac.FE << "\t" << formfac.FM;
      g_formfacs_fe.at( j )->SetPoint( g_formfacs_fe.at( j )->GetN(), q2, formfac.FE );
      g_formfacs_fm.at( j )->SetPoint( g_formfacs_fm.at( j )->GetN(), q2, formfac.FM );
      ++j;
    }
    out << "\n";
  }
  CG_INFO( "main" )
    << "Scan written in \"" << output_file << "\".";
  out.close();

  vector<int> colours = { kBlack, kRed+1, kBlue-2, kGreen+1, kOrange+1 };

  for ( auto& plt : map<const char*,vector<TGraph*> >{ { "FE", g_formfacs_fe }, { "FM", g_formfacs_fm } } ) {
    cepgen::Canvas c( plt.first, Form( "M_{X} = %g GeV/c^{2}", mx ) );
    c.SetLogy();
    TMultiGraph mg;
    size_t i = 0;
    for ( auto& gr : plt.second ) {
      mg.Add( gr );
      gr->SetLineColor( colours.at( i ) );
      c.AddLegendEntry( gr, gr->GetTitle(), "l" );
      ++i;
    }
    mg.Draw( "al" );
    mg.GetHistogram()->SetTitle( Form( ";Q^{2};%s", plt.first ) );
    c.Prettify( mg.GetHistogram() );
    c.Save( "pdf" );
  }

  return 0;
}
