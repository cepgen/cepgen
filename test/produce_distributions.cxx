#include "CepGen/Generator.h"
#include "CepGen/Cards/ConfigReader.h"

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
  const unsigned long num_gen_events = 1e4;

  if ( argc<2 ) {
    InError( Form( "Usage: %s [input card]", argv[0] ) );
    return -1;
  }
  mg.setParameters( CepGen::Cards::ConfigReader( argv[1] ).parameters() );

  TH1D h_mass( "invm", "Dilepton invariant mass\\d#sigma/dm\\GeV?.2f", 1000, 0., 500. ),
       h_ptpair( "ptpair", "Dilepton p_{T}\\d#sigma/dp_{T}\\GeV?.2f", 500, 0., 50. );

  std::ostringstream gen_name;
  gen_name << mg.parameters->process()->name();
  Information( Form( "Process name: %s", gen_name.str().c_str() ) );
  //mg.parameters->taming_functions.dump();

  for ( unsigned int i=0; i<num_gen_events; i++ ) {
    CepGen::Event* ev = mg.generateOneEvent();
    if ( i%100==0 ) Information( Form( "Produced event #%d", i ) );
    const auto pl1 = ev->getOneByRole( CepGen::Particle::CentralParticle1 ).momentum(),
               pl2 = ev->getOneByRole( CepGen::Particle::CentralParticle2 ).momentum();
    h_mass.Fill( ( pl1+pl2 ).mass() );
    h_ptpair.Fill( ( pl1+pl2 ).pt() );
  }
  const double weight = mg.crossSection()/num_gen_events;
  h_mass.Scale( weight );
  h_ptpair.Scale( weight );

  produce_plot( "dilepton_invm", &h_mass );
  produce_plot( "dilepton_ptpair", &h_ptpair );

  return 0;
}
