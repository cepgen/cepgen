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
  const float min_xbj = 1.e-3, max_xbj = 1.0, q2 = ( argc>1 ) ? atof( argv[1] ) : 2.5;
  const char* q2_str = ( argc>2 ) ? argv[2] : "2.5";
  const unsigned int npoints = 1000;

  TGraph g_sy_fe_100, g_sy_fm_100;
  TGraph g_fb_fe_100, g_fb_fm_100;
  TGraph g_su_fe_100, g_su_fm_100;

  const float mp2 = pow( CepGen::Particle::massFromPDGId( CepGen::Particle::Proton ), 2 );

  for ( unsigned int i=0; i<npoints; i++ ) {
    const float xbj = min_xbj + i*( max_xbj-min_xbj )/( npoints-1 );
    const float mx2 = mp2 + ( 1.-xbj )*q2;

    CepGen::FormFactors ff_sy = CepGen::FormFactors::SuriYennie( q2, mp2, mx2 );
    g_sy_fe_100.SetPoint( i, xbj, ff_sy.FE );
    g_sy_fm_100.SetPoint( i, xbj, ff_sy.FM );

    CepGen::FormFactors ff_fb = CepGen::FormFactors::FioreBrasse( q2, mp2, mx2 );
    g_fb_fe_100.SetPoint( i, xbj, ff_fb.FE );
    g_fb_fm_100.SetPoint( i, xbj, ff_fb.FM );

    CepGen::FormFactors ff_su = CepGen::FormFactors::SzczurekUleshchenko( q2, mp2, mx2 );
    g_su_fe_100.SetPoint( i, xbj, ff_su.FE );
    g_su_fm_100.SetPoint( i, xbj, ff_su.FM );
  }

  CepGen::Canvas c( "test", Form( "CepGen proton form factors, Q^{2} = %s GeV^{2}", q2_str ) );
  c.SetLegendX1( 0.4 );

  TMultiGraph mg;

  g_sy_fe_100.SetLineWidth( 3 );
  mg.Add( &g_sy_fe_100, "l" );
  c.AddLegendEntry( &g_sy_fe_100, "Suri-Yennie, F_{E}", "l" );

  g_sy_fm_100.SetLineStyle( 2 );
  g_sy_fm_100.SetLineWidth( 3 );
  mg.Add( &g_sy_fm_100, "l" );
  c.AddLegendEntry( &g_sy_fm_100, "Suri-Yennie, F_{M}", "l" );

  g_fb_fe_100.SetLineColor( kRed+1 );
  g_fb_fe_100.SetLineWidth( 3 );
  mg.Add( &g_fb_fe_100, "l" );
  c.AddLegendEntry( &g_fb_fe_100, "Fiore-Brasse, F_{E}", "l" );

  g_fb_fm_100.SetLineStyle( 2 );
  g_fb_fm_100.SetLineColor( kRed+1 );
  g_fb_fm_100.SetLineWidth( 3 );
  mg.Add( &g_fb_fm_100, "l" );
  c.AddLegendEntry( &g_fb_fm_100, "Fiore-Brasse, F_{M}", "l" );

  g_su_fe_100.SetLineColor( kGreen+2 );
  g_su_fe_100.SetLineWidth( 3 );
  mg.Add( &g_su_fe_100, "l" );
  c.AddLegendEntry( &g_su_fe_100, "Szczurek-Uleshchenko, F_{E}", "l" );

  g_su_fm_100.SetLineStyle( 2 );
  g_su_fm_100.SetLineColor( kGreen+2 );
  g_su_fm_100.SetLineWidth( 3 );
  mg.Add( &g_su_fm_100, "l" );
  c.AddLegendEntry( &g_su_fm_100, "Szczurek-Uleshchenko, F_{M}", "l" );

  mg.Draw( "alpr" );
  mg.SetTitle( "x_{Bj}\\Proton form factor" );

  c.Prettify( mg.GetHistogram() );
  //mg.GetYaxis()->SetRangeUser( 1.e-10, 10. );
  //mg.GetYaxis()->SetRangeUser( 0., 0.6 );
  mg.GetXaxis()->SetLimits( min_xbj, max_xbj );
  //c.SetLogx();
  //c.SetLogy();

  c.Save( "pdf" );

  return 0;
}
