#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "Canvas.h"

#include <iostream>

#include "CepGen/Utils/Functional.h"
#include "CepGen/Utils/ArgumentsParser.h"

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

  TGraph gr_fb, gr_rt;
  TF1 f_rt( "f_rt", "TMath::Min(1.,TMath::Exp(-x/10))", min_x, max_x );

  cepgen::utils::Functional<1> test1( "min(1,exp(-x/10))", { { "x" } } );
  for ( unsigned short i = 0; i < num_points; ++i ) {
    const double x = min_x + ( max_x-min_x )/( num_points-1 )*i;
    gr_rt.SetPoint( i, x, f_rt.Eval( x ) );
    try {
      gr_fb.SetPoint( i, x, test1.eval( x ) );
    } catch ( const cepgen::Exception& e ) { e.dump(); }
  }
  TGraph gr_diff;
  double chi2 = 0.;
  for ( unsigned short i = 0; i < gr_fb.GetN(); ++i ) {
    gr_diff.SetPoint( i, gr_fb.GetX()[i], gr_fb.GetY()[i]-gr_rt.GetY()[i] );
    chi2 += pow( gr_fb.GetY()[i]-gr_rt.GetY()[i], 2 );
  }
  chi2 = sqrt( chi2 );
  if ( chi2 > 1.e-9 )
    throw CG_FATAL( "main" ) << "Test failed with chi2 = " << chi2 << "!";
  cout << "Test passed!" << endl;

  if ( draw ) {
    cepgen::Canvas c( "test_graph", "CepGen validation", true );
    TMultiGraph mg;
    mg.Add( &gr_fb );
    mg.Add( &gr_rt );
    mg.Draw( "al" );
    c.AddLegendEntry( &gr_fb, "Functional", "l" );
    c.AddLegendEntry( &gr_rt, "ROOT", "l" );
    gr_fb.SetLineWidth( 3 );
    gr_fb.SetLineStyle( 2 );
    c.Prettify( mg.GetHistogram() );
    auto ratio = c.RatioPlot( gr_fb.GetHistogram(), gr_rt.GetHistogram(), -1., 1., "l" );
    ratio->SetLineColor( kRed );
    ratio->SetLineWidth( 3 );
    gr_diff.SetLineStyle( 2 );
    gr_diff.SetLineColor( kBlue );
    gr_diff.Draw( "same" );
    c.Save( "pdf" );
  }
  return 0;
}
