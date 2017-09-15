#include "CepGen/Export/HepMCHandler.h"

using namespace std;
using namespace CepGen;

int main() {

  OutputHandler::HepMCHandler writer( "example.dat" );
  writer.setCrossSection(1., 2.);

  Event ev;
  
  Particle p1( Particle::IncomingBeam1, Particle::Proton );
  p1.setMomentum( 1., -15., 100. );
  p1.setStatus( Particle::Incoming );
  ev.addParticle(p1);

  Particle p2( Particle::IncomingBeam2, Particle::Electron );
  p2.setMomentum( 10., 5., 3200. );
  p2.setStatus( Particle::Incoming );
  ev.addParticle(p2);

  ev.dump();
  
  writer << &ev;

  return 0;
}
