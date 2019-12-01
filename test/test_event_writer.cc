#include "CepGen/Generator.h"

#include "CepGen/Physics/PDG.h"
#include "CepGen/Event/Event.h"

#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Utils/ArgumentsParser.h"

#include <iostream>

using namespace std;
using namespace cepgen;

int main( int argc, char* argv[] )
{
  Generator gen;

  string type;
  bool list;

  ArgumentsParser( argc, argv )
    .addOptionalArgument( "format", "type of format to build", "hepmc", &type )
    .addOptionalArgument( "list", "list all formats", false, &list, 'l' )
    .parse();

  if ( list ) {
    cout
      << "List of export modules available:\n"
      << "=================================\n";
    for ( const auto& mod : io::ExportModuleFactory::get().modules() )
      cout << mod << std::endl;
    return 0;
  }

  auto writer = io::ExportModuleFactory::get().build( type );
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
