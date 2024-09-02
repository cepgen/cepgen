/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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
#include "CepGenAddOns/Common/EventUtils.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::initialise();

  string type;
  bool list;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("format,f", "type of format to build", &type, "hepmc")
      .addOptionalArgument("list,l", "list all formats", &list, false)
      .parse();

  if (list) {
    CG_LOG.log([](auto& log) {
      log << "List of export modules available:\n"
          << "=================================";
      for (const auto& mod : cepgen::EventExporterFactory::get().modules())
        log << "\n" << mod;
    });
    return 0;
  }

  auto writer = cepgen::EventExporterFactory::get().build(type);
  writer->setCrossSection(cepgen::Value{1., 2.});

  {  // first test: simple event content
    cepgen::Event ev;

    cepgen::Particle p1(cepgen::Particle::Role::IncomingBeam1, cepgen::PDG::proton);
    p1.setMomentum(1., -15., 100.);
    p1.setStatus(cepgen::Particle::Status::Incoming);
    ev.addParticle(p1);

    cepgen::Particle p2(cepgen::Particle::Role::IncomingBeam2, cepgen::PDG::electron);
    p2.setMomentum(10., 5., 3200.);
    p2.setStatus(cepgen::Particle::Status::Incoming);
    ev.addParticle(p2);

    const auto ev_old = ev;
    CG_DEBUG("main") << "Event content:\n" << ev;
    *writer << ev;
    CG_TEST(ev == ev_old, "[" + type + "] event content preservation by output (simple event)");
  }
  {  // second test: simple event content with parentage
    auto ev = cepgen::Event::minimal();
    const auto ev_old = ev;
    CG_DEBUG("main") << "Event content:\n" << ev;
    *writer << ev;
    CG_TEST(ev == ev_old, "[" + type + "] event content preservation by output (simple+parentage event)");
  }
  {  // third test: realistic lpair event content
    auto ev = cepgen::utils::generateLPAIREvent();
    const auto ev_old = ev;
    CG_DEBUG("main") << "Event content:\n" << ev;
    *writer << ev;
    CG_TEST(ev == ev_old, "[" + type + "] event content preservation by output (lpair event)");
  }

  CG_TEST_SUMMARY;
}
