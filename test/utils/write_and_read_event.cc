/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/EventUtils.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::Generator gen;

  const auto writers = cepgen::EventExporterFactory::get().modules(),
             readers = cepgen::EventImporterFactory::get().modules();
  vector<string> common;
  for (const auto& mod : writers)
    if (cepgen::utils::contains(readers, mod))
      common.emplace_back(mod);

  cepgen::ArgumentsParser(argc, argv).addOptionalArgument("modules,m", "modules to test", &common, common).parse();
  CG_INFO("main") << "Will test with the following writer/reader pairs: " << common << ".";

  const auto evt_base = cepgen::utils::generateLPAIREvent();
  const auto cross_section = cepgen::Value{42.4242, 0.4242};

  CG_DEBUG("main") << "Input event to be tested:\n" << evt_base;

  for (const auto& mod : common) {
    string temp_file = "output.txt";  // default value

    {  // write event to output file
      auto writer = cepgen::EventExporterFactory::get().build(mod);
      temp_file = writer->parameters().get<string>("filename");
      writer->initialise(gen.runParameters());
      writer->setCrossSection(cross_section);
      const auto wrote = ((*writer) << evt_base);
      CG_TEST(wrote, "event export: " + mod);
      if (!wrote)
        continue;
    }
    {  // read back output file
      auto reader =
          cepgen::EventImporterFactory::get().build(mod, cepgen::ParametersList().set<string>("filename", temp_file));
      reader->initialise(gen.runParameters());
      cepgen::Event evt_in;
      CG_TEST_EQUAL(((*reader) >> evt_in), true, "event re-import: " + mod);
      CG_TEST_EQUAL(evt_in.size(), evt_base.size(), "event re-import size: " + mod);
      CG_TEST_EQUAL(reader->crossSection(), cross_section, "stored cross-section: " + mod);
      for (const auto& role : {cepgen::Particle::Role::IncomingBeam1,
                               cepgen::Particle::Role::IncomingBeam2,
                               cepgen::Particle::Role::OutgoingBeam1,
                               cepgen::Particle::Role::OutgoingBeam2,
                               cepgen::Particle::Role::Parton1,
                               cepgen::Particle::Role::Parton2}) {
        ostringstream os_role;
        os_role << role;
        CG_TEST_EQUAL(evt_in.oneWithRole(role).integerPdgId(),
                      evt_base.oneWithRole(role).integerPdgId(),
                      "PDG of " + os_role.str() + ": " + mod);
        CG_TEST_EQUIV(evt_in.oneWithRole(role).momentum().px(),
                      evt_base.oneWithRole(role).momentum().px(),
                      "x-momentum of " + os_role.str() + ": " + mod);
        CG_TEST_EQUIV(evt_in.oneWithRole(role).momentum().py(),
                      evt_base.oneWithRole(role).momentum().py(),
                      "y-momentum of " + os_role.str() + ": " + mod);
        CG_TEST_EQUIV(evt_in.oneWithRole(role).momentum().pz(),
                      evt_base.oneWithRole(role).momentum().pz(),
                      "z-momentum of " + os_role.str() + ": " + mod);
      }
    }
  }
  CG_TEST_SUMMARY;
}
