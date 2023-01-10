/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace cepgen;

int main(int argc, char* argv[]) {
  initialise();

  string type;
  bool list;

  ArgumentsParser(argc, argv)
      .addOptionalArgument("format", "type of format to build", &type, "hepmc")
      .addOptionalArgument("list,l", "list all formats", &list, false)
      .parse();

  if (list) {
    CG_LOG.log([](auto& log) {
      log << "List of export modules available:\n"
          << "=================================";
      for (const auto& mod : EventExporterFactory::get().modules())
        log << "\n" << mod;
    });
    return 0;
  }

  auto writer = EventExporterFactory::get().build(type);
  writer->setCrossSection(1., 2.);

  Event ev;

  Particle p1(Particle::IncomingBeam1, PDG::proton);
  p1.setMomentum(1., -15., 100.);
  p1.setStatus(Particle::Status::Incoming);
  ev.addParticle(p1);

  Particle p2(Particle::IncomingBeam2, PDG::electron);
  p2.setMomentum(10., 5., 3200.);
  p2.setStatus(Particle::Status::Incoming);
  ev.addParticle(p2);

  ev.dump();

  *writer << ev;

  CG_TEST_SUMMARY;
}
