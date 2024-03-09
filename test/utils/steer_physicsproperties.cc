#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Generator.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/Test.h"

class TestObject : public cepgen::SteeredObject<TestObject> {
public:
  explicit TestObject(const cepgen::ParametersList& params = cepgen::ParametersList())
      : cepgen::SteeredObject<TestObject>(params),
        particle_props_(steer<cepgen::ParticleProperties>("particleProps")) {}
  const cepgen::ParticleProperties& particleProperties() const { return particle_props_; }
  static cepgen::ParametersDescription description() {
    auto desc = cepgen::ParametersDescription();
    desc.addAs<int, cepgen::pdgid_t>("particleProps", cepgen::PDG::muon);
    return desc;
  }

private:
  cepgen::ParticleProperties particle_props_;
};

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  {
    //const auto pdesc = cepgen::ParametersList().setAs<int, cepgen::pdgid_t>("particleProps", cepgen::PDG::muon);
    auto pprop = cepgen::ParticleProperties{};
    pprop.fermion = true;
    pprop.name = "laurenton";
    pprop.pdgid = 42;
    pprop.mass = 42.4242;
    const auto pdesc = cepgen::ParametersList().set<cepgen::ParticleProperties>("particleProps", pprop);

    auto object = TestObject(pdesc);

    CG_TEST_EQUAL(cepgen::PDG::get()(42), pprop, "Full particle properties object");
  }
  {
    cepgen::initialise();
    TestObject object;
    CG_TEST_EQUAL(
        object.particleProperties(), cepgen::PDG::get()(cepgen::PDG::muon), "Full particle properties object");
  }
  CG_TEST_SUMMARY;
}
