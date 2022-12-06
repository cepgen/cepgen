#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/Message.h"

int main() {
  cepgen::initialise();

  auto paw = cepgen::io::ExportModuleFactory::get().build("paw");
  paw->initialise(cepgen::Parameters());

  auto evt = cepgen::Event();
  auto& p1 = evt.addParticle(cepgen::Particle::Role::CentralSystem, 22).get();
  p1.setMomentum(cepgen::Momentum(1, 2, 3, 4));

  (*paw) << evt;

  return 0;
}
