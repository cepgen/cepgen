#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"

using namespace std;

int main() {
  cepgen::initialise();

  auto pythia = cepgen::EventModifierFactory::get().build("pythia6");

  cepgen::Event evt;
  auto tau = cepgen::Particle(cepgen::Particle::CentralSystem, 15, cepgen::Particle::Status::Undecayed);
  tau.setMomentum(0., 0., 1000.);
  evt.addParticle(tau);

  double weight;
  pythia->run(evt, weight, true);

  CG_LOG << evt;

  return 0;
}
