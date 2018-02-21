#include "CepGen/Generator.h"
#include "CepGen/Cards/PythiaHandler.h"

#include "Canvas.h"
#include "TH1.h"

#include <sstream>

void produce_plot( const char* name, TH1* hist )
{
  CepGen::Canvas c( name, "CepGen Simulation" );
  hist->Draw( "hist" );
  c.Prettify( hist );
  c.SetLogy();
  c.Save( "pdf" );
}

int main( int argc, char* argv[] )
{
  CepGen::Generator mg;
  //CepGen::Logger::get().level = CepGen::Logger::Debug;

  if ( argc < 2 ) {
    InError( Form( "Usage: %s [input card]", argv[0] ) );
    return -1;
  }
  mg.setParameters( CepGen::Cards::PythiaHandler( argv[1] ).parameters() );

  TH1D h_mass( "invm", "Dilepton invariant mass\\d#sigma/dM\\GeV?.2f", 1000, 0., 500. ),
       h_ptpair( "ptpair", "Dilepton p_{T}\\d#sigma/dp_{T}\\GeV?.2f", 500, 0., 50. ),
       h_ptsingle( "pt_single", "Single lepton p_{T}\\d#sigma/dp_{T}\\?.2f", 100, 0., 100. ),
       h_etasingle( "eta_single", "Single lepton #eta\\d#sigma/d#eta\\?.2f", 60, -3., 3. );

  std::ostringstream gen_name;
  gen_name << mg.parameters->process()->name();
  Information( Form( "Process name: %s", gen_name.str().c_str() ) );
  //mg.parameters->taming_functions.dump();

  for ( unsigned int i = 0; i < mg.parameters->generation.maxgen; ++i ) {
    const CepGen::Event& ev = *mg.generateOneEvent();
    if ( i%100==0 ) Information( Form( "Produced event #%d", i ) );
    const auto central_system = ev.getByRole( CepGen::Particle::CentralSystem );
    const auto pl1 = central_system[0].momentum(), pl2 = central_system[1].momentum();
    h_mass.Fill( ( pl1+pl2 ).mass() );
    h_ptpair.Fill( ( pl1+pl2 ).pt() );
    h_ptsingle.Fill( pl1.pt() );
    h_etasingle.Fill( pl1.eta() );
  }
  const double weight = mg.crossSection()/mg.parameters->generation.maxgen;
  h_mass.Scale( weight );
  h_ptpair.Scale( weight );

  produce_plot( "dilepton_invm", &h_mass );
  produce_plot( "dilepton_ptpair", &h_ptpair );
  produce_plot( "singlelepton_pt", &h_ptsingle );
  produce_plot( "singlelepton_eta", &h_etasingle );

  return 0;
}
