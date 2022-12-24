#include <HepMC/GenEvent.h>
#include <HepMC/IO_GenEvent.h>

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main() {
  cepgen::initialise();

  auto evt = cepgen::Event();

  auto p1 = cepgen::Particle(cepgen::Particle::Role::CentralSystem, 2212, cepgen::Particle::Status::FinalState);
  p1.momentum().setP(1., 2., 3., 4.);
  evt.addParticle(p1);

  auto p2 = cepgen::Particle(cepgen::Particle::Role::CentralSystem, 2212, cepgen::Particle::Status::FinalState);
  p2.momentum().setP(2., 4., 6., 8.);
  evt.addParticle(p2);

  auto temp_file = "/tmp/test_hepmc.out";
  {
    auto hepmc_out = cepgen::EventExporterFactory::get().build(
        "hepmc2", cepgen::ParametersList().set<string>("filename", temp_file));
    (*hepmc_out) << evt;
  }
  {
    auto hepmc_in = cepgen::EventImporterFactory::get().build("hepmc2");
    HepMC::IO_GenEvent reader(temp_file, std::ios::in);
    HepMC::GenEvent event;
    CG_TEST_EQUAL(reader.fill_next_event(&event), true, "Event re-import [HepMC2]");

    auto evt_in = hepmc_in->convert(event);
    CG_TEST_EQUAL(evt_in.size(), evt.size(), "Event re-import size");
    for (const auto& part : evt_in.particles()) {
      //CG_TEST_EQUAL(part.pdgId(), evt[part.id()].pdgId(), "Event re-import");
    }
  }

  CG_TEST_SUMMARY;
}
