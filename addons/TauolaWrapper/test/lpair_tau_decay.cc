/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventModifier.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/EventUtils.h"
#include "CepGen/Utils/Test.h"

using namespace std;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  string rng_name;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("random-generator,r", "type of random number generator to use", &rng_name, ""s)
      .parse();
  cepgen::initialise();

  auto tauola_parameters = cepgen::ParametersList();
  if (!rng_name.empty())  // a particular random number generator is specified
    tauola_parameters.set("randomGenerator", cepgen::ParametersList().setName(rng_name));

  auto tauola = cepgen::EventModifierFactory::get().build("tauola", tauola_parameters);
  if (!tauola) {
    CG_LOG << "Failed to retrieve the Tauola interface!";
    return -1;
  }
  tauola->initialise(cepgen::RunParameters());

  auto event = cepgen::utils::generateLPAIREvent();
  // modify two-lepton system kinematics to generate taus
  auto oc = event[cepgen::Particle::Role::CentralSystem];
  oc[0].get().setPdgId(cepgen::PDG::tau, -1);
  oc[0].get().momentum().setMass(cepgen::PDG::get().mass(cepgen::PDG::tau));
  oc[1].get().setPdgId(cepgen::PDG::tau, +1);
  oc[1].get().momentum().setMass(cepgen::PDG::get().mass(cepgen::PDG::tau));

  double weight = 1.;
  const auto event_size_before_decay = event.size();
  tauola->run(event, weight);
  CG_LOG << event;

  CG_TEST(event.size() != event_size_before_decay, "decay was performed");

  CG_TEST_SUMMARY;
}
