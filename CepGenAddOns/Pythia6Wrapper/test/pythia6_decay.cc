#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/Hadroniser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main() {
  cepgen::initialise();

  auto pythia = cepgen::EventModifierFactory::get().build("pythia6");

  cepgen::Event evt;
  auto tau = cepgen::Particle(cepgen::Particle::CentralSystem, 15, cepgen::Particle::Status::Undecayed);
  tau.setMomentum(cepgen::Momentum(0., 0., 1000.), false);
  evt.addParticle(tau);
  const auto evt_size_bef = evt.size();

  double weight;
  pythia->run(evt, weight, false);

  CG_LOG << evt;
  CG_TEST_EQUAL(evt[0].status(), cepgen::Particle::Status::Resonance, "tau 'decayed' status");
  CG_TEST(evt_size_bef != evt.size(), "decay");

  CG_TEST_SUMMARY;
}
