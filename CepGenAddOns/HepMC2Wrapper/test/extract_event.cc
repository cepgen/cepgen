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

  auto evt = cepgen::Event();
  {
    auto p0 =
        cepgen::Particle(cepgen::Particle::Role::IncomingBeam1, 2212, cepgen::Particle::Status::PrimordialIncoming);
    p0.momentum().setP(0.5, 1., 1.5, 2.);
    evt.addParticle(p0);
    auto p1 =
        cepgen::Particle(cepgen::Particle::Role::IncomingBeam2, 2212, cepgen::Particle::Status::PrimordialIncoming);
    p1.momentum().setP(1., 2., 3., 4.);
    evt.addParticle(p1);
    auto p2 = cepgen::Particle(cepgen::Particle::Role::OutgoingBeam1, 2212, cepgen::Particle::Status::FinalState);
    p2.momentum().setP(2., 4., 6., 8.);
    evt.addParticle(p2);
    p2.addMother(p0);
    auto p3 = cepgen::Particle(cepgen::Particle::Role::OutgoingBeam2, 2212, cepgen::Particle::Status::FinalState);
    p3.momentum().setP(4., 8., 12., 16.);
    evt.addParticle(p3);
    p3.addMother(p1);
    auto p4 = cepgen::Particle(cepgen::Particle::Role::CentralSystem, 12, cepgen::Particle::Status::FinalState);
    p4.momentum().setP(8., 16., 24., 32.);
    evt.addParticle(p4);
    p4.addMother(p0);
    p4.addMother(p1);
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
    for (const auto& part : evt_in.particles()) {
      //CG_TEST_EQUAL(part.pdgId(), evt[part.id()].pdgId(), "Event re-import");
    }
  }

  CG_TEST_SUMMARY;
}
