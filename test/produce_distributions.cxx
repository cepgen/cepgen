#include "CepGen/Generator.h"
#include "CepGen/Cards/Handler.h"

#include "Canvas.h"
#include "TH1.h"

#include <sstream>

void produce_plot( const char* name, TH1* hist )
{
  CepGen::Canvas c( name, "CepGen Simulation" );
  hist->Draw();
  c.Prettify( hist );
  c.Save( "pdf" );
}

int main( int argc, char* argv[] )
{
  CepGen::Generator mg;

  if ( argc<2 ) {
    InError( Form( "Usage: %s [input card]", argv[0] ) );
    return -1;
  }
  mg.setParameters( CepGen::Cards::LpairReader( argv[1] ).parameters() );

  TH1D h_mass( "invm", "Dilepton invariant mass\\Events\\GeV?.2f", 1000, 0., 100. ),
       h_ptpair( "ptpair", "Dilepton p_{T}\\Events\\GeV?.2f", 500, 0., 50. );

  std::ostringstream gen_name;
  gen_name << mg.parameters->process->name();
  Information( Form( "Process name: %s", gen_name.str().c_str() ) );

  for ( unsigned int i=0; i<1e4; i++ ) {
    CepGen::Event* ev = mg.generateOneEvent();
    if ( i%100==0 ) Information( Form( "Produced event #%d", i ) );
    const CepGen::Particle::Momentum pl1 = ev->getOneByRole( CepGen::Particle::CentralParticle1 )->momentum(),
                                     pl2 = ev->getOneByRole( CepGen::Particle::CentralParticle2 )->momentum();
    h_mass.Fill( ( pl1+pl2 ).mass() );
    h_ptpair.Fill( ( pl1+pl2 ).pt() );
  }

  produce_plot( "dilepton_invm", &h_mass );
  produce_plot( "dilepton_ptpair", &h_ptpair );

  return 0;
}
