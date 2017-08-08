#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/Particle.h"
#include "Canvas.h"

#include "TGraph.h"
#include "TMultiGraph.h"

#include <iostream>

using namespace std;

int
main( int argc, char* argv[] )
{
  const float min_q2 = 1., max_q2 = 1.e5, mx2 = 1.e4;
  const unsigned int npoints = 1000;

  TGraph g_sy_fe_100, g_sy_fm_100;
  TGraph g_fb_fe_100, g_fb_fm_100;
  TGraph g_su_fe_100, g_su_fm_100;

  const float mp2 = pow( Particle::GetMassFromPDGId( Particle::Proton ), 2 );

  for ( unsigned int i=0; i<npoints; i++ ) {
    const float q2 = min_q2 + i*( max_q2-min_q2 )/(npoints-1);

    FormFactors ff_sy = SuriYennieFormFactors( q2, mp2, mx2 );
    g_sy_fe_100.SetPoint( i, q2, ff_sy.FE );
    g_sy_fm_100.SetPoint( i, q2, ff_sy.FM );

    FormFactors ff_fb = FioreBrasseFormFactors( q2, mp2, mx2 );
    g_fb_fe_100.SetPoint( i, q2, ff_fb.FE );
    g_fb_fm_100.SetPoint( i, q2, ff_fb.FM );

    FormFactors ff_su = SzczurekUleschenkoFormFactors( q2, mp2, mx2 );
    g_su_fe_100.SetPoint( i, q2, ff_su.FE );
    g_su_fm_100.SetPoint( i, q2, ff_su.FM );
  }

  Canvas c( "test" );
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

  g_su_fe_100.SetLineColor( kGreen-1 );
  g_su_fe_100.SetLineWidth( 3 );
  mg.Add( &g_su_fe_100, "l" );
  c.AddLegendEntry( &g_su_fe_100, "Szczurek-Uleschenko, F_{E}", "l" );

  g_su_fm_100.SetLineStyle( 2 );
  g_su_fm_100.SetLineColor( kGreen-1 );
  g_su_fm_100.SetLineWidth( 3 );
  mg.Add( &g_su_fm_100, "l" );
  c.AddLegendEntry( &g_su_fm_100, "Szczurek-Uleschenko, F_{M}", "l" );

  mg.Draw( "alpr" );
  mg.SetTitle( "Q^{2}\\Proton form factor" );

  c.Prettify( mg.GetHistogram() );
  mg.GetYaxis()->SetRangeUser( 1.e-10, 10. );
  mg.GetXaxis()->SetLimits( min_q2, max_q2 );
  c.SetLogy();

  c.Save( "pdf" );

  return 0;
}
