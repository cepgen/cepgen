#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/Particle.h"
#include "test/Canvas.h"

#include "TGraph.h"
#include "TMultiGraph.h"

#include <iostream>

using namespace std;

int
main( int argc, char* argv[] )
{
  const float min_xbj = 1.e-5, max_xbj = 0.99, q2 = ( argc>1 ) ? atof( argv[1] ) : 2.5;
  const char* q2_str = ( argc>2 ) ? argv[2] : std::to_string( q2 ).c_str();
  const unsigned int npoints = 5000;

  TGraph g_sy_f1, g_sy_f2;
  TGraph g_fb_f1, g_fb_f2;
  TGraph g_su_f1, g_su_f2;

  const bool use_logarithmic_x = ( argc>3 ) ? atoi( argv[3] ) : false;

  for ( unsigned int i=0; i<npoints; i++ ) {
    float xbj;
    if ( use_logarithmic_x ) {
      const float min_lxbj = log10( min_xbj ), max_lxbj = log10( max_xbj );
      xbj = pow( 10, min_lxbj + i*( max_lxbj-min_lxbj )/( npoints-1 ) );
      std::cout << min_lxbj << "\t" << max_lxbj << "\t" << xbj << std::endl;
    }
    else xbj = min_xbj + i*( max_xbj-min_xbj )/( npoints-1 );

    auto sf_sy = CepGen::StructureFunctions::SuriYennie( q2, xbj ),
         sf_fb = CepGen::StructureFunctions::FioreBrasse( q2, xbj ),
         sf_su = CepGen::StructureFunctions::SzczurekUleshchenko( q2, xbj );

    g_sy_f1.SetPoint( i, xbj, sf_sy.F1 );
    g_sy_f2.SetPoint( i, xbj, sf_sy.F2 );

    g_fb_f1.SetPoint( i, xbj, sf_fb.F1 );
    g_fb_f2.SetPoint( i, xbj, sf_fb.F2 );

    g_su_f1.SetPoint( i, xbj, sf_su.F1 );
    g_su_f2.SetPoint( i, xbj, sf_su.F2 );
    std::cout << sf_fb << std::endl;
  }

  CepGen::Canvas c( "test", Form( "CepGen proton structure functions, Q^{2} = %s GeV^{2}", q2_str ) );
  c.SetLegendX1( 0.4 );

  TMultiGraph mg;

  /*g_sy_f1.SetLineWidth( 3 );
  mg.Add( &g_sy_f1, "l" );
  c.AddLegendEntry( &g_sy_f1, "Suri-Yennie, F_{1}", "l" );*/

  g_sy_f2.SetLineStyle( 2 );
  g_sy_f2.SetLineWidth( 3 );
  mg.Add( &g_sy_f2, "l" );
  c.AddLegendEntry( &g_sy_f2, "Suri-Yennie, F_{2}", "l" );

  /*g_fb_f1.SetLineColor( kRed+1 );
  g_fb_f1.SetLineWidth( 3 );
  mg.Add( &g_fb_f1, "l" );
  c.AddLegendEntry( &g_fb_f1, "Fiore-Brasse, F_{1}", "l" );*/

  //g_fb_f2.SetLineStyle( 2 );
  g_fb_f2.SetLineColor( kRed+1 );
  g_fb_f2.SetLineWidth( 3 );
  mg.Add( &g_fb_f2, "l" );
  c.AddLegendEntry( &g_fb_f2, "Fiore-Brasse, F_{2}", "l" );

  /*g_su_f1.SetLineColor( kGreen+2 );
  g_su_f1.SetLineWidth( 3 );
  mg.Add( &g_su_f1, "l" );
  c.AddLegendEntry( &g_su_f1, "Szczurek-Uleshchenko, F_{1}", "l" );*/

  //g_su_f2.SetLineStyle( 2 );
  g_su_f2.SetLineColor( kGreen+2 );
  g_su_f2.SetLineWidth( 3 );
  mg.Add( &g_su_f2, "l" );
  c.AddLegendEntry( &g_su_f2, "Szczurek-Uleshchenko, F_{2}", "l" );

  mg.Draw( "alpr" );
  mg.SetTitle( "x_{Bj}\\Proton form factor" );

  c.Prettify( mg.GetHistogram() );
  //mg.GetYaxis()->SetRangeUser( 1.e-10, 10. );
  //mg.GetYaxis()->SetRangeUser( 0., 0.6 );
  mg.GetXaxis()->SetLimits( min_xbj, max_xbj );
  if ( use_logarithmic_x ) c.SetLogx();
  //c.SetLogy();

  c.Save( "pdf" );

  return 0;
}
