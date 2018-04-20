#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "Canvas.h"

#include <iostream>

#include "CepGen/Core/Functional.h"

using namespace std;

int main( int argc, char* argv[] )
{
  const unsigned short num_points = 100;
  const double min_x = -1., max_x = 10.;
  TGraph gr_fb, gr_rt;
  TF1 f_rt( "f_rt", "TMath::Min(1.,TMath::Exp(-x/10))", min_x, max_x );

  CepGen::Functional<1> test1( "min(1,exp(-x/10))", { { "x" } } );
  for ( unsigned short i = 0; i < num_points; ++i ) {
    const double x = min_x + ( max_x-min_x )/( num_points-1 )*i;
    gr_rt.SetPoint( i, x, f_rt.Eval( x ) );
    try {
      gr_fb.SetPoint( i, x, test1.eval( x ) );
    } catch ( const CepGen::Exception& e ) { e.dump(); }
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

  if ( argc > 1 && !strcmp( argv[1], "draw" ) ) {
    CepGen::Canvas c( "test_graph", "CepGen validation" );
    TMultiGraph mg;
    mg.Add( &gr_fb );
    mg.Add( &gr_rt );
    mg.Add( &gr_diff );
    mg.Draw( "al" );
    c.AddLegendEntry( &gr_fb, "Functional", "l" );
    c.AddLegendEntry( &gr_rt, "ROOT", "l" );
    c.AddLegendEntry( &gr_diff, "Difference", "l" );
    gr_fb.SetLineWidth( 3 );
    gr_fb.SetLineStyle( 2 );
    gr_diff.SetLineColor( kRed );
    c.Prettify( mg.GetHistogram() );
    c.Save( "pdf" );
  }
  return 0;
}
