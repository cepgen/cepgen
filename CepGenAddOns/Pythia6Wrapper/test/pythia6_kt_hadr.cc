#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"

using namespace std;

cepgen::Event generate_event() {
  auto evt = cepgen::Event::minimal(2);
  auto& ib1 = evt.oneWithRole(cepgen::Particle::IncomingBeam1);
  ib1.setPdgId(cepgen::PDG::proton);
  ib1.setMomentum(cepgen::Momentum::fromPxPyPzE(0., 0., 6.5e3, -1.), false);
  auto& ib2 = evt.oneWithRole(cepgen::Particle::IncomingBeam2);
  ib2.setPdgId(cepgen::PDG::proton);
  ib2.setMomentum(cepgen::Momentum::fromPxPyPzE(0., 0., -6.5e3, -1.), false);
  auto& ob1 = evt.oneWithRole(cepgen::Particle::OutgoingBeam1);
  ob1.setPdgId(cepgen::PDG::proton);
  ob1.setMomentum(cepgen::Momentum::fromPxPyPzE(-7.875321, 8.186351, 6.403512e3, 6.403704e3), true);
  auto& ob2 = evt.oneWithRole(cepgen::Particle::OutgoingBeam2);
  ob2.setPdgId(cepgen::PDG::proton);
  ob2.setMomentum(cepgen::Momentum::fromPxPyPzE(-2.725610e-2, 7.565269e-3, -6.425336e3, 6.425336e3), false);
  auto& parton1 = evt.oneWithRole(cepgen::Particle::Parton1);
  parton1.setPdgId(cepgen::PDG::photon);
  parton1.setMomentum(cepgen::Momentum::fromPxPyPzE(7.875321, -8.186351, 9.648800e1, 9.629600e1), true);
  auto& parton2 = evt.oneWithRole(cepgen::Particle::Parton2);
  parton2.setPdgId(cepgen::PDG::photon);
  parton2.setMomentum(cepgen::Momentum::fromPxPyPzE(2.725610e-2, -7.565269e-3, -7.466409e1, 7.466409e1), true);
  evt.oneWithRole(cepgen::Particle::Intermediate).setMomentum(parton1.momentum() + parton2.momentum(), true);
  auto oc = evt[cepgen::Particle::CentralSystem];
  oc[0].get().setPdgId(cepgen::PDG::muon);
  oc[0].get().setMomentum(cepgen::Momentum::fromPxPyPzE(2.193109e1, -6.725967e1, -4.248568e1, 8.252200e1), false);
  oc[1].get().setPdgId(cepgen::PDG::muon);
  oc[1].get().setMomentum(cepgen::Momentum::fromPxPyPzE(-1.402852e1, 5.906575e1, 6.430959e1, 8.843809e1), false);
  return evt;
}

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  auto gen = cepgen::Generator();

  auto evt = generate_event();
  evt.dump();

  gen.parametersRef().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set(
          "kinematics",
          cepgen::ParametersList()
              .set<double>("cmEnergy", 13.e3)
              .setAs<int, cepgen::mode::Kinematics>("mode", cepgen::mode::Kinematics::InelasticElastic))));

  auto cg_pythia = cepgen::EventModifierFactory::get().build("pythia6");
  cg_pythia->setCrossSection(cepgen::Value{1.46161e-1, 1.25691e-3});
  cg_pythia->initialise(gen.parametersRef());
  double evt_weight = 1.;

  const auto evt_before_particles = evt.particles().size();
  cg_pythia->run(evt, evt_weight, true);
  CG_TEST(evt_weight == 1., "no event weight modification in fast mode");
  CG_TEST(evt.particles().size() == evt_before_particles, "no event modification in fast mode");

  cg_pythia->run(evt, evt_weight, false);

  CG_DEBUG("main") << "Pythia 6-filtered event:\n" << evt;

  CG_TEST_EQUAL(evt_weight, 1., "event weight");
  CG_TEST(evt(cepgen::Particle::Role::OutgoingBeam1).size() > 1, "decayed diffractive beam system");
  CG_TEST(evt(cepgen::Particle::Role::OutgoingBeam2).size() == 1, "undecayed elastic beam system");
  cepgen::Momentum daugh_total_momentum;
  for (const auto& daugh : evt.stableDaughters(evt(cepgen::Particle::Role::OutgoingBeam1)[0], true))
    daugh_total_momentum += daugh.get().momentum();
  CG_TEST_EQUIV((daugh_total_momentum - evt(cepgen::Particle::Role::OutgoingBeam1)[0].momentum()).p(),
                0.,
                "diffractive system momentum balance");

  CG_TEST_SUMMARY;
}
