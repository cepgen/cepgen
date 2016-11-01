#include "core/MCGen.h"

#include "Canvas.h"
#include "TH1.h"

#include <sstream>

void produce_plot( const char* name, TH1* hist )
{
  Canvas c( name, "CepGen Simulation" );
  hist->Draw();
  c.Prettify( hist );
  c.Save( "pdf" );
}

int main( int argc, char* argv[] )
{
  MCGen mg;

  if ( argc<2 ) {
    InError( Form( "Usage: %s [input card]", argv[0] ) );
    return -1;
  }
  if ( !mg.parameters->ReadConfigFile( argv[1] ) ) {
    InError( Form( "Error reading the configuration!\n\t"
                   "Please check your input file (%s)", argv[1] ) );
    return -1;
  }

  TH1D h_mass( "invm", "Dilepton invariant mass\\Events\\GeV?.2f", 1000, 0., 100. ),
       h_ptpair( "ptpair", "Dilepton p_{T}\\Events\\GeV?.2f", 500, 0., 50. );

  std::ostringstream gen_name;
  gen_name << mg.parameters->process;
  Information( Form( "Process name: %s", gen_name.str().c_str() ) );

  for ( unsigned int i=0; i<1e5; i++ ) {
    Event* ev = mg.GenerateOneEvent();
    if ( i%100==0 ) Information( Form( "Produced event #%d", i ) );
    const Particle::Momentum pl1 = ev->GetOneByRole( Particle::CentralParticle1 )->GetMomentum(),
                             pl2 = ev->GetOneByRole( Particle::CentralParticle2 )->GetMomentum();
    h_mass.Fill( ( pl1+pl2 ).M() );
    h_ptpair.Fill( ( pl1+pl2 ).Pt() );
  }

  produce_plot( "dilepton_invm", &h_mass );
  produce_plot( "dilepton_ptpair", &h_ptpair );

  return 0;
}
