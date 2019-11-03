#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Core/Exception.h"

#include "Canvas.h"
#include "TH1.h"

#include "ArgumentsParser.h"

using namespace std;

unique_ptr<TH1D> h_mass, h_ptpair, h_ptsingle, h_etasingle;
void process_event( const cepgen::Event& ev, unsigned long event_id )
{
  cout << event_id << endl;
  const auto& central_system = ev[cepgen::Particle::CentralSystem];
  const auto& pl1 = central_system[0].momentum(), pl2 = central_system[1].momentum();
  h_mass->Fill( ( pl1+pl2 ).mass() );
  h_ptpair->Fill( ( pl1+pl2 ).pt() );
  h_ptsingle->Fill( pl1.pt() );
  h_etasingle->Fill( pl1.eta() );
}

int main( int argc, char* argv[] )
{
  cepgen::Generator mg;

  string card;
  int num_events;

  cepgen::ArgumentsParser( argc, argv )
    .addArgument( "input", "input card", &card, 'i' )
    .addOptionalArgument( "num-events", "number of events to generate", -1, &num_events, 'n' )
    .parse();

  mg.setParameters( cepgen::card::Handler::parse( card.c_str() ) );
  if ( num_events >= 0 ) { // user specified a number of events to generate
    mg.parameters().generation().maxgen = num_events;
    mg.parameters().generation().enabled = num_events > 0;
  }

  h_mass.reset( new TH1D( "invm", ";Dilepton invariant mass;d#sigma/dM (pb/GeV)", 500, 0., 500. ) );
  h_ptpair.reset( new TH1D( "ptpair", ";Dilepton p_{T};d#sigma/dp_{T} (pb/GeV)", 500, 0., 50. ) );
  h_ptsingle.reset( new TH1D( "pt_single", ";Single lepton p_{T};d#sigma/dp_{T} (pb/GeV)", 100, 0., 100. ) );
  h_etasingle.reset( new TH1D( "eta_single", ";Single lepton #eta;d#sigma/d#eta (pb)\\?.2f", 60, -3., 3. ) );

  CG_INFO( "main" ) << "Process name: " << mg.parameters().processName() << ".";
  //mg.parameters->taming_functions.dump();

  mg.generate( process_event );

  const double weight = mg.crossSection()/mg.parameters().generation().maxgen;
  h_mass->Scale( weight, "width" );
  h_ptpair->Scale( weight, "width" );
  h_ptsingle->Scale( weight, "width" );
  h_etasingle->Scale( weight, "width" );

  const unordered_map<const char*, TH1D&> plots = {
    { "dilepton_invm", *h_mass },
    { "dilepton_ptpair", *h_ptpair },
    { "singlelepton_pt", *h_ptsingle },
    { "singlelepton_eta", *h_etasingle },
  };
  for ( auto& plt : plots ) {
    cepgen::Canvas c( plt.first, "CepGen Simulation" );
    plt.second.Draw( "hist" );
    c.Prettify( &plt.second );
    c.SetLogy();
    c.Save( "pdf" );
  }

  return 0;
}
