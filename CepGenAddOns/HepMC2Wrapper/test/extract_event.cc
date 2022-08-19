#include <HepMC/GenEvent.h>

#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/Test.h"

int main() {
  cepgen::initialise();

  auto evt = cepgen::Event();

  auto p1 = cepgen::Particle(cepgen::Particle::Role::CentralSystem, 2212, cepgen::Particle::Status::FinalState);
  p1.momentum().setP(1., 2., 3., 4.);
  evt.addParticle(p1);

  auto p2 = cepgen::Particle(cepgen::Particle::Role::CentralSystem, 2212, cepgen::Particle::Status::FinalState);
  p2.momentum().setP(2., 4., 6., 8.);
  evt.addParticle(p2);

  auto hepmc_out = cepgen::EventExporterFactory::get().build("hepmc");
  (*hepmc_out) << evt;

  auto hepmc_in = cepgen::EventImporterFactory::get().build("hepmc");
  HepMC::GenEvent event;
  event.print();
  hepmc_in->convert(&event, evt);
  CG_LOG << evt;

  CG_TEST_SUMMARY;
}
