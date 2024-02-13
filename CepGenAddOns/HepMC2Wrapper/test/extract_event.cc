#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  auto evt = cepgen::Event::minimal(2);
  {
    auto &ip1 = evt.oneWithRole(cepgen::Particle::Role::IncomingBeam1),
         &ip2 = evt.oneWithRole(cepgen::Particle::Role::IncomingBeam2),
         &op1 = evt.oneWithRole(cepgen::Particle::Role::OutgoingBeam1),
         &op2 = evt.oneWithRole(cepgen::Particle::Role::OutgoingBeam2),
         &part1 = evt.oneWithRole(cepgen::Particle::Role::Parton1),
         &part2 = evt.oneWithRole(cepgen::Particle::Role::Parton2);
    ip1.setPdgId((long)2212);
    ip1.momentum().setP(0.5, 1., 1.5, 2.);
    ip2.setPdgId((long)2212);
    ip2.momentum().setP(1., 2., 3., 4.);
    op1.setPdgId((long)2212);
    op1.momentum().setP(2., 4., 6., 8.);
    op2.setPdgId((long)2212);
    op2.momentum().setP(4., 8., 12., 16.);
    part1.setPdgId((long)22);
    part2.setPdgId((long)22);
    evt[cepgen::Particle::Role::CentralSystem][0].get().momentum().setP(8., 16., 24., 32.);
    evt[cepgen::Particle::Role::CentralSystem][1].get().momentum().setP(16., 32., 64., 128.);
  }

  auto temp_file = "/tmp/test_hepmc.out";
  {
    auto hepmc_out = cepgen::EventExporterFactory::get().build(
        "hepmc2", cepgen::ParametersList().set<string>("filename", temp_file));
    (*hepmc_out) << evt;
  }
  {
    auto hepmc_in = cepgen::EventImporterFactory::get().build(
        "hepmc2", cepgen::ParametersList().set<string>("filename", temp_file));
    cepgen::Event evt_in;
    CG_TEST_EQUAL(((*hepmc_in) >> evt_in), true, "Event re-import [HepMC2]");
    CG_TEST_EQUAL(evt_in.size(), evt.size(), "Event re-import size");
    for (const auto& role : {cepgen::Particle::Role::IncomingBeam1,
                             cepgen::Particle::Role::IncomingBeam2,
                             cepgen::Particle::Role::OutgoingBeam1,
                             cepgen::Particle::Role::OutgoingBeam2,
                             cepgen::Particle::Role::Parton1,
                             cepgen::Particle::Role::Parton2}) {
      ostringstream os_role;
      os_role << role;
      CG_TEST_EQUAL(evt_in.oneWithRole(role).pdgId(), evt.oneWithRole(role).pdgId(), "PDG of " + os_role.str());
      CG_TEST_EQUAL(evt_in.oneWithRole(role).momentum().px(),
                    evt.oneWithRole(role).momentum().px(),
                    "x-momentum of " + os_role.str());
      CG_TEST_EQUAL(evt_in.oneWithRole(role).momentum().py(),
                    evt.oneWithRole(role).momentum().py(),
                    "y-momentum of " + os_role.str());
      CG_TEST_EQUAL(evt_in.oneWithRole(role).momentum().pz(),
                    evt.oneWithRole(role).momentum().pz(),
                    "z-momentum of " + os_role.str());
    }
  }

  CG_TEST_SUMMARY;
}
