#include "CepGen/Core/ExportHandler.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Event/Event.h"

using namespace std;
using namespace cepgen;

int main()
{
  auto writer = io::ExportHandler::get().build( "hepmc" );
  writer->setCrossSection( 1., 2. );

  Event ev;

  Particle p1( Particle::IncomingBeam1, PDG::proton );
  p1.setMomentum( 1., -15., 100. );
  p1.setStatus( Particle::Status::Incoming );
  ev.addParticle(p1);

  Particle p2( Particle::IncomingBeam2, PDG::electron );
  p2.setMomentum( 10., 5., 3200. );
  p2.setStatus( Particle::Status::Incoming );
  ev.addParticle(p2);

  ev.dump();

  *writer << ev;

  return 0;
}
