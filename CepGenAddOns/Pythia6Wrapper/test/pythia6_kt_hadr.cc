#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/Common/EventUtils.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  auto gen = cepgen::Generator();

  auto evt = cepgen::utils::generateLPAIREvent();
  evt.dump();

  gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set(
          "kinematics",
          cepgen::ParametersList()
              .set<double>("cmEnergy", 13.e3)
              .setAs<int, cepgen::mode::Kinematics>("mode", cepgen::mode::Kinematics::InelasticElastic))));

  auto cg_pythia = cepgen::EventModifierFactory::get().build("pythia6");
  cg_pythia->setCrossSection(cepgen::Value{1.46161e-1, 1.25691e-3});
  cg_pythia->initialise(gen.runParameters());
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
