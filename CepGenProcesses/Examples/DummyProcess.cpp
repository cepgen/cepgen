/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"

using namespace cepgen;

auto make_pdgids_pair = [](pdgid_t pair) {
  // from one PDG identifier, return a pair of particle/anti-particle
  return spdgids_t{(spdgid_t)pair, -(spdgid_t)pair};
};

/// Compute a dummy factorised matrix element
class DummyProcess final : public cepgen::proc::FactorisedProcess {
public:
  explicit DummyProcess(const cepgen::ParametersList& params)
      : FactorisedProcess(params, make_pdgids_pair(params.get<ParticleProperties>("pair").pdgid)),
        value_(steer<double>("value")),  // example on how to steer a parameter from user configuration
        flag_(steer<int>("flag")) {}

  static ParametersDescription description() {
    // This static member allows to generate a human- and machine-readable description of this algorithm/process ;
    // for instance, switches and flags can be documented using one of the DocumentationGenerator objects to
    // generate e.g. HTML, or text descriptions.
    auto desc = cepgen::proc::FactorisedProcess::description();
    desc.setDescription("Dummy 2-to-4 process");
    desc.add<double>("value", 1.).setDescription("a floating point value given by the user");
    desc.add<int>("flag", 1)
        .setDescription(
            "another value given by the user, e.g. to switch between several modes that can be handled by this "
            "fantastic example")
        .allow(1, "value is equal to 1")
        .allow(42, "value is equal to 42");
    return desc;
  }

private:
  void prepareFactorisedPhaseSpace() override {
    // This method allows you to prepare the matrix element computation with the kinematics information
    // retrieved from the `kin_` member inherited from the cepgen::proc::Process base object.
    const auto& cs_prop = PDG::get()(psgen_->central().at(0));
    CG_DEBUG("DummyProcess:prepare") << "Produced particles: " << cs_prop.descr << " ("
                                     << "mass = " << cs_prop.mass << " GeV.";
  }
  double computeFactorisedMatrixElement() override {
    // This is where the central, parton-factorised matrix element is being computed.
    return value_;
  }

  const double value_;
  const int flag_;
};

REGISTER_PROCESS("dummy", DummyProcess);  // register the process into the runtime database
