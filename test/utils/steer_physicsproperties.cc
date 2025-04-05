#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Generator.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/Test.h"

using namespace std::string_literals;

class TestObject : public cepgen::SteeredObject<TestObject> {
public:
  explicit TestObject(const cepgen::ParametersList& params = cepgen::ParametersList())
      : cepgen::SteeredObject<TestObject>(params),
        particle_props_(steer<cepgen::ParticleProperties>("particleProps")) {}
  const cepgen::ParticleProperties& particleProperties() const { return particle_props_; }
  static cepgen::ParametersDescription description() {
    auto desc = cepgen::ParametersDescription();
    desc.addAs<cepgen::pdgid_t>("particleProps", cepgen::PDG::muon);
    return desc;
  }

private:
  const cepgen::ParticleProperties particle_props_;
};

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::initialise();
  {
    auto particle_properties = cepgen::ParticleProperties{};
    particle_properties.fermion = true;
    particle_properties.name = "laurenton"s;
    particle_properties.pdgid = 42;
    particle_properties.mass = 42.4242;
    particle_properties.charges = {-3, 3};

    auto object =
        TestObject(cepgen::ParametersList().set<cepgen::ParticleProperties>("particleProps", particle_properties));
    CG_TEST_EQUAL(cepgen::PDG::get()(42), particle_properties, "part.prop. registered in PDG database");
    CG_TEST_EQUAL(object.particleProperties(), particle_properties, "part.prop. retrieved from steered object");
  }
  {
    TestObject object;
    CG_TEST_EQUAL(object.particleProperties().pdgid, cepgen::PDG::muon, "part. default from steered object");
    CG_TEST_EQUAL(object.particleProperties().mass,
                  cepgen::PDG::get()(cepgen::PDG::muon).mass,
                  "part.prop. default from steered object");
  }
  CG_TEST_SUMMARY;
}
