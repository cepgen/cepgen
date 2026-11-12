#include <CepGen/Core/RunParameters.h>
#include <CepGen/Event/Event.h>
#include <CepGen/EventFilter/EventModifier.h>
#include <CepGen/Generator.h>
#include <CepGen/Modules/EventModifierFactory.h>
#include <CepGen/Modules/ProcessFactory.h>
#include <CepGen/Physics/Modes.h>
#include <CepGen/Process/Process.h>
#include <CepGen/Utils/ArgumentsParser.h>
#include <CepGen/Utils/EventUtils.h>
#include <CepGen/Utils/Test.h>

#include <cmath>

using namespace std::string_literals;

int main(int argc, char* argv[]) {
  int num_events;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("num-events,n", "number of events to hadronise", &num_events, 100)
      .parse();

  auto gen = cepgen::Generator();
  gen.runParameters().setProcess(cepgen::ProcessFactory::get().build(
      "lpair",
      cepgen::ParametersList().set(
          "kinematics",
          cepgen::ParametersList()
              .set<double>("cmEnergy", 13.e3)
              .setAs<int, cepgen::mode::Kinematics>("mode", cepgen::mode::Kinematics::InelasticElastic))));

  auto hadroniser = cepgen::EventModifierFactory::get().build("sherpa");
  hadroniser->initialise(gen.runParameters());

  // events recycled by the generator: system ancestry should not leak from one event to another
  auto evt = cepgen::utils::generateLPAIREvent();
  evt.freeze();
  const auto system_role = cepgen::Particle::Role::OutgoingBeam1;
  const auto system_momentum = evt(system_role)[0].momentum();
  bool all_balanced = true, children_consistent = true, weights_unchanged = true, all_hadronised = true,
       all_fragmented = true;
  for (int i = 0; i < num_events; ++i) {
    evt.restore();
    double weight = -1.;
    all_hadronised &= hadroniser->run(evt, weight, false);
    weights_unchanged &= weight == 1.;
    cepgen::Momentum total;
    for (const auto& stable : evt.stableChildren(evt(system_role)[0], true))
      total += stable.get().momentum();
    all_balanced &= (total - system_momentum).p() < 1.e-6;
    const auto& system = evt(system_role)[0];
    all_fragmented &= system.status() == cepgen::Particle::Status::Fragmented;
    for (const auto& child : system.children())
      children_consistent &= evt(child).mothers().count(system.id()) == 1;
  }
  CG_TEST(all_hadronised, "all systems hadronised");
  CG_TEST(weights_unchanged, "unit event weights");
  CG_TEST(all_fragmented, "system flagged as fragmented");
  CG_TEST(all_balanced, "momentum balance in recycled events");
  CG_TEST(children_consistent, "no stale filiation in recycled events");

  evt.restore();
  evt[evt(system_role)[0].id()].setMomentum(cepgen::Momentum::fromPxPyPzM(0., 0., 100., 1.05), true);
  double weight = -1.;
  CG_TEST(!hadroniser->run(evt, weight, true), "veto of a system below threshold in fast mode");

  CG_TEST_SUMMARY;
}
