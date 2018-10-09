#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Generator.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Core/Exception.h"

#include "Canvas.h"
#include "TH1.h"

#include <sstream>

using namespace std;

void produce_plot( const char* name, TH1* hist )
{
  cepgen::Canvas c( name, "CepGen Simulation" );
  hist->Draw( "hist" );
  c.Prettify( hist );
  c.SetLogy();
  c.Save( "pdf" );
}

unique_ptr<TH1D> h_mass, h_ptpair, h_ptsingle, h_etasingle;
void process_event( const cepgen::Event& ev, unsigned long event_id )
{
  const auto central_system = ev.getByRole( cepgen::Particle::CentralSystem );
  const auto pl1 = central_system[0].momentum(), pl2 = central_system[1].momentum();
  h_mass->Fill( ( pl1+pl2 ).mass() );
  h_ptpair->Fill( ( pl1+pl2 ).pt() );
  h_ptsingle->Fill( pl1.pt() );
  h_etasingle->Fill( pl1.eta() );
}

int main( int argc, char* argv[] )
{
  cepgen::Generator mg;

  if ( argc < 2 )
    throw CG_FATAL( "main" ) << "Usage: " << argv[0] << " [input card]";
  mg.setParameters( cepgen::card::PythonHandler( argv[1] ).parameters() );

  h_mass = unique_ptr<TH1D>( new TH1D( "invm", ";Dilepton invariant mass;d#sigma/dM (pb/GeV)", 500, 0., 500. ) );
  h_ptpair = unique_ptr<TH1D>( new TH1D( "ptpair", ";Dilepton p_{T};d#sigma/dp_{T} (pb/GeV)", 500, 0., 50. ) );
  h_ptsingle = unique_ptr<TH1D>( new TH1D( "pt_single", ";Single lepton p_{T};d#sigma/dp_{T} (pb/GeV)", 100, 0., 100. ) );
  h_etasingle = unique_ptr<TH1D>( new TH1D( "eta_single", ";Single lepton #eta;d#sigma/d#eta (pb)\\?.2f", 60, -3., 3. ) );

  CG_INFO( "main" ) << "Process name: " << mg.parameters->processName() << ".";
  //mg.parameters->taming_functions.dump();

  mg.generate( process_event );

  const double weight = mg.crossSection()/mg.parameters->generation.maxgen;
  h_mass->Scale( weight, "width" );
  h_ptpair->Scale( weight, "width" );
  h_ptsingle->Scale( weight, "width" );
  h_etasingle->Scale( weight, "width" );

  produce_plot( "dilepton_invm", h_mass.get() );
  produce_plot( "dilepton_ptpair", h_ptpair.get() );
  produce_plot( "singlelepton_pt", h_ptsingle.get() );
  produce_plot( "singlelepton_eta", h_etasingle.get() );

  return 0;
}
