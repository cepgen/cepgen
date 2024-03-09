#include "CepGen/Core/Exception.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Test.h"
#include "CepGenAddOns/MadGraphWrapper/MadGraphProcess.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  const cepgen::pdgid_t my_part = 13;
  const auto my_part_mass = 42.;

  cepgen::PDG::get().define(cepgen::ParticleProperties(my_part, "la", "laurentino", 0., my_part_mass, 0., 3, true));
  CG_LOG << cepgen::PDG::get()(my_part);

  auto mg5 = cepgen::ProcessFactory::get().build(
      "mg5_aMC",
      cepgen::ParametersList()
          .set("kinematicsGenerator", cepgen::ParametersList().setName("coll2to4"s))
          .set("extraParticles", cepgen::ParametersList().set("la", cepgen::PDG::get()(my_part)))
          .set("process", "a a > la+ la-"s));
  mg5->initialise();

  const auto& proc_evt = mg5->event();
  CG_TEST_EQUAL(proc_evt.oneWithRole(cepgen::Particle::Role::Parton1).pdgId(), cepgen::PDG::photon, "parton 1 PDG id");
  CG_TEST_EQUAL(proc_evt.oneWithRole(cepgen::Particle::Role::Parton2).pdgId(), cepgen::PDG::photon, "parton 2 PDG id");
  CG_TEST_EQUAL(proc_evt(cepgen::Particle::Role::CentralSystem).size(), 2, "cent.part.multiplicity");
  for (size_t i = 0; i < 2; ++i) {
    const auto& cent = proc_evt(cepgen::Particle::Role::CentralSystem).at(i);
    CG_TEST_EQUAL(cent.integerPdgId(), short((i == 0 ? -1 : +1) * my_part), "cent." + to_string(i) + " PDG id");
    CG_TEST_EQUAL(cent.momentum().mass(), my_part_mass, "cent." + to_string(i) + " mass");
  }

  CG_TEST_SUMMARY;
}
