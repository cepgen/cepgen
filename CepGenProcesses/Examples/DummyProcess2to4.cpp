/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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
#include "CepGen/Process/Process2to4.h"

using namespace cepgen;

/// Compute a dummy 2-to-4 matrix element
class DummyProcess2to4 final : public cepgen::proc::Process2to4 {
public:
  explicit DummyProcess2to4(const cepgen::ParametersList& params)
      : Process2to4(params, {PDG::photon, PDG::photon}, steer<ParticleProperties>("pair").pdgid),
        value_(steer<double>("value")) {}

  static ParametersDescription description() {
    auto desc = cepgen::proc::Process2to4::description();
    desc.setDescription("Dummy 2-to-4 process (kt-factor.)");
    desc.add<double>("value", 1.);
    return desc;
  }

private:
  void prepareProcessKinematics() override {
    // this method allows you to prepare the matrix element computation with the kinematics information
    // retrieved from the `kin_` member inherited from the cepgen::proc::Process base object.
    const auto& cs_prop = PDG::get()(produced_parts_.at(0));
    CG_DEBUG("DummyProcess2to4:prepare") << "Produced particles: " << cs_prop_.descr << " ("
                                         << "mass = " << cs_prop_.mass << " GeV.";
  }
  double computeCentralMatrixElement() const override { return value_; }

  const double value_;
};

// register process
REGISTER_PROCESS("dummy", DummyProcess2to4);
