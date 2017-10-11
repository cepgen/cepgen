#include "CepGen/Event/Particle.h"
#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
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

  TGraph g_sy_f2, g_fb_f2, g_su_f2, g_bdh_f2, g_cteq_f2, g_mrst_f2;
  TGraph g_allm97_f2, g_allm_hht_f2, g_allm_hht_ft_f2;
  TGraph g_lux_f2;
  TGraph g_cb_f2;

  const bool use_logarithmic_x = ( argc>3 ) ? atoi( argv[3] ) : false;

  /*CepGen::SF::GenericLHAPDF cteq( "cteq6l1" );
  CepGen::SF::GenericLHAPDF mrst( "MRST2004qed_proton" );
  CepGen::SF::GenericLHAPDF lux( "LUXqed17_plus_PDF4LHC15_nnlo_100" );*/

  for ( unsigned int i=0; i<npoints; i++ ) {
    float xbj;
    if ( use_logarithmic_x ) {
      const float min_lxbj = log10( min_xbj ), max_lxbj = log10( max_xbj );
      xbj = pow( 10, min_lxbj + i*( max_lxbj-min_lxbj )/( npoints-1 ) );
    }
    else xbj = min_xbj + i*( max_xbj-min_xbj )/( npoints-1 );

    auto sf_sy = CepGen::StructureFunctionsBuilder::get( CepGen::StructureFunctions::SuriYennie, q2, xbj ),
         sf_fb = CepGen::StructureFunctionsBuilder::get( CepGen::StructureFunctions::FioreBrasse, q2, xbj ),
         sf_su = CepGen::StructureFunctionsBuilder::get( CepGen::StructureFunctions::SzczurekUleshchenko, q2, xbj ),
         sf_allm97 = CepGen::StructureFunctionsBuilder::get( CepGen::StructureFunctions::ALLM97, q2, xbj ),
         //sf_allm_hht = CepGen::SF::ALLM( q2, xbj, CepGen::SF::ALLMParameterisation::hht_allm() ),
         //sf_allm_hht_ft = CepGen::SF::ALLM( q2, xbj, CepGen::SF::ALLMParameterisation::hht_allm_ft() ),
         //sf_bdh = CepGen::SF::BlockDurandHa( q2, xbj ),
         //sf_cteq = cteq( q2, xbj ),
         //sf_mrst = mrst( q2, xbj ),
         //sf_lux = lux( q2, xbj ),
         sf_cb = CepGen::StructureFunctionsBuilder::get( CepGen::StructureFunctions::ChristyBosted, q2, xbj );

    g_sy_f2.SetPoint( i, xbj, sf_sy.F2 );
    g_fb_f2.SetPoint( i, xbj, sf_fb.F2 );
    g_su_f2.SetPoint( i, xbj, sf_su.F2 );
    //g_bdh_f2.SetPoint( i, xbj, sf_bdh.F2 );
    //g_cteq_f2.SetPoint( i, xbj, sf_cteq.F2 );
    //g_mrst_f2.SetPoint( i, xbj, sf_mrst.F2 );
    //g_lux_f2.SetPoint( i, xbj, sf_lux.F2 );
    g_cb_f2.SetPoint( i, xbj, sf_cb.F2 );

    g_allm97_f2.SetPoint( i, xbj, sf_allm97.F2 );
    //g_allm_hht_f2.SetPoint( i, xbj, sf_allm_hht.F2 );
    //g_allm_hht_ft_f2.SetPoint( i, xbj, sf_allm_hht_ft.F2 );
  }

  CepGen::Canvas c( "test", Form( "CepGen proton structure functions, Q^{2} = %s GeV^{2}", q2_str ) );
  c.SetLegendX1( 0.4 );

  TMultiGraph mg;

  /*g_sy_f2.SetLineStyle( 2 );
  g_sy_f2.SetLineWidth( 3 );
  mg.Add( &g_sy_f2, "l" );
  c.AddLegendEntry( &g_sy_f2, "Suri-Yennie", "l" );*/

  //g_fb_f2.SetLineStyle( 2 );
  g_fb_f2.SetLineColor( kRed+1 );
  g_fb_f2.SetLineWidth( 3 );
  mg.Add( &g_fb_f2, "l" );
  c.AddLegendEntry( &g_fb_f2, "Fiore-Brasse", "l" );

  //g_su_f2.SetLineStyle( 2 );
  /*g_su_f2.SetLineColor( kGreen+2 );
  g_su_f2.SetLineWidth( 3 );
  mg.Add( &g_su_f2, "l" );
  c.AddLegendEntry( &g_su_f2, "Szczurek-Uleshchenko", "l" );*/

  g_allm97_f2.SetLineColor( kBlue+1 );
  g_allm97_f2.SetLineWidth( 3 );
  mg.Add( &g_allm97_f2, "l" );
  c.AddLegendEntry( &g_allm97_f2, "Abramowicz et al. 97", "l" );

  /*g_allm_hht_f2.SetLineColor( kBlue+1 );
  g_allm_hht_f2.SetLineWidth( 3 );
  g_allm_hht_f2.SetLineStyle( 2 );
  mg.Add( &g_allm_hht_f2, "l" );
  c.AddLegendEntry( &g_allm_hht_f2, "Abramowicz et al. HHT", "l" );

  g_allm_hht_ft_f2.SetLineColor( kBlue+1 );
  g_allm_hht_ft_f2.SetLineWidth( 3 );
  g_allm_hht_ft_f2.SetLineStyle( 3 );
  mg.Add( &g_allm_hht_ft_f2, "l" );
  c.AddLegendEntry( &g_allm_hht_ft_f2, "Abramowicz et al. HHT-FT", "l" );*/

  g_cb_f2.SetLineColor( kMagenta );
  g_cb_f2.SetLineWidth( 3 );
  mg.Add( &g_cb_f2, "l" );
  c.AddLegendEntry( &g_cb_f2, "Christy-Bosted", "l" );

  /*g_bdh_f2.SetLineColor( kOrange );
  g_bdh_f2.SetLineWidth( 3 );
  mg.Add( &g_bdh_f2, "l" );
  c.AddLegendEntry( &g_bdh_f2, "Block-Durand-Ha", "l" );

  g_cteq_f2.SetLineColor( kMagenta );
  g_cteq_f2.SetLineWidth( 3 );
  mg.Add( &g_cteq_f2, "l" );
  c.AddLegendEntry( &g_cteq_f2, "CTEQ6", "l" );

  g_mrst_f2.SetLineColor( kMagenta+1 );
  g_mrst_f2.SetLineWidth( 3 );
  mg.Add( &g_mrst_f2, "l" );
  c.AddLegendEntry( &g_mrst_f2, "MRST2004 (QED/proton)", "l" );*/

  /*g_lux_f2.SetLineColor( kOrange );
  g_lux_f2.SetLineWidth( 3 );
  mg.Add( &g_lux_f2, "l" );
  c.AddLegendEntry( &g_lux_f2, "LUXqed", "l" );*/

  mg.Draw( "alpr" );
  mg.SetTitle( "x_{Bj}\\Proton form factor F_{2}" );

  c.Prettify( mg.GetHistogram() );
  //mg.GetYaxis()->SetRangeUser( 1.e-10, 10. );
  mg.GetYaxis()->SetRangeUser( 0., 0.8 );
  mg.GetXaxis()->SetLimits( min_xbj, max_xbj );
  if ( use_logarithmic_x ) c.SetLogx();
  //c.SetLogy();

  c.Save( "pdf" );

  return 0;
}
