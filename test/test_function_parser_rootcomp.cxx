#include "CepGen/Modules/FunctionalFactory.h"
#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Core/Exception.h"

#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "Canvas.h"

#include <iostream>

using namespace std;

int main( int argc, char* argv[] )
{
  bool draw;
  int num_points;
  double min_x, max_x;

  cepgen::ArgumentsParser( argc, argv )
    .addOptionalArgument( "draw", "do draw the canvas?", false, &draw, 'd' )
    .addOptionalArgument( "num-points", "number of points to consider", 100, &num_points, 'n' )
    .addOptionalArgument( "min-x", "minimal range", -1., &min_x, 'l' )
    .addOptionalArgument( "max-x", "maximal range", 1., &max_x, 'H' )
    .parse();

  TGraph gr_rt;
  TF1 f_rt( "f_rt", "TMath::Min(1.,TMath::Exp(-x/10))", min_x, max_x );
  for ( unsigned short i = 0; i < num_points; ++i ) {
    const double x = min_x + ( max_x-min_x )/( num_points-1 )*i;
    gr_rt.SetPoint( i, x, f_rt.Eval( x ) );
  }

  map<string, TGraph> m_gr_fb, m_gr_diff;
  for ( const auto& func : cepgen::utils::FunctionalFactory::get().modules() ) {
    CG_LOG( "main" ) << "Building \"" << func << "\" functional.";
    try {
      auto test = cepgen::utils::FunctionalFactory::get().build( func, cepgen::ParametersList()
        .set<std::string>( "expression", "min(1,exp(-x/10))" )
        .set<std::vector<std::string> >( "variables", { "x" } ) );
      double chi2 = 0;
      for ( unsigned short i = 0; i < num_points; ++i ) {
        const double x = min_x + ( max_x-min_x )/( num_points-1 )*i;
        const double val = (*test)( x ), val_ref = gr_rt.GetY()[i];
        m_gr_fb[func].SetPoint( i, x, val );
        m_gr_diff[func].SetPoint( i, x, val-val_ref );
        chi2 += pow( val-val_ref, 2 );
      }
      chi2 = sqrt( chi2 );
      if ( chi2 > 1.e-9 )
        throw CG_FATAL( "main" ) << "Test failed with chi2 = " << chi2 << "!";
    } catch ( const cepgen::Exception& ) {
      CG_WARNING( "main" )
        << "Exception encountered in \"" << func << "\" functional builder.";
      continue;
    }
  }

  cout << "Test passed!" << endl;

  if ( draw ) {
    cepgen::Canvas c( "test_graph", "CepGen validation", true );
    TMultiGraph mg;
    mg.Add( &gr_rt );
    c.AddLegendEntry( &gr_rt, "ROOT", "l" );
    size_t i = 0;
    for ( auto& gr_fb : m_gr_fb ) {
      mg.Add( &gr_fb.second );
      gr_fb.second.SetLineWidth( 3 );
      gr_fb.second.SetLineStyle( 2+( i++ ) );
      c.AddLegendEntry( &gr_fb.second, Form( "Functional (%s)", gr_fb.first.c_str() ), "l" );
    }
    i = 0;
    for ( auto& gr_diff : m_gr_diff ) {
      gr_diff.second.SetLineStyle( 2+( i++ ) );
      gr_diff.second.SetLineColor( kBlue );
      gr_diff.second.Draw( "same" );
    }
    mg.Draw( "al" );
    c.Prettify( mg.GetHistogram() );
    c.Save( "pdf" );
    i = 0;
    for ( auto& gr_fb : m_gr_fb ) {
      auto ratio = c.RatioPlot( gr_fb.second.GetHistogram(), gr_rt.GetHistogram(), -1., 1., ( i == 0 ? "al" : "l,same" ) );
      ratio->SetLineColor( kRed );
      ratio->SetLineWidth( 3 );
      ratio->SetLineStyle( 2+( i++ ) );
    }
  }
  return 0;
}
