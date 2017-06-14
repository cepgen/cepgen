#include "export/HepMCHandler.h"

using namespace std;

int main() {

  OutputHandler::HepMCHandler writer( "example.dat" );
  writer.setCrossSection(1., 2.);

  Event ev;
  
  Particle p1( Particle::IncomingBeam1, Particle::Proton );
  p1.SetMomentum( 1., -15., 100. );
  p1.status = Particle::Incoming;
  ev.AddParticle(p1);

  Particle p2( Particle::IncomingBeam2, Particle::Electron );
  p2.SetMomentum( 10., 5., 3200. );
  p2.status = Particle::Incoming;
  ev.AddParticle(p2);

  ev.Dump();
  
  writer << &ev;

  return 0;
}
