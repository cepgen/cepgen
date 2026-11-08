/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

  const auto tau1_momentum = cepgen::Momentum::fromPxPyPzM(0., 0., +100., cepgen::PDG::get().mass(cepgen::PDG::tau)),
             tau2_momentum = cepgen::Momentum::fromPxPyPzM(0., 0., -100., cepgen::PDG::get().mass(cepgen::PDG::tau));

  cepgen::Event ev;
  cepgen::Particle pho(cepgen::Particle::Role::CentralSystem, cepgen::PDG::photon, cepgen::Particle::Status::Resonance);
  pho.setMomentum(tau1_momentum + tau2_momentum);
  ev.addParticle(pho);
  cepgen::Particle tau1(
      cepgen::Particle::Role::CentralSystem, +(cepgen::spdgid_t)cepgen::PDG::tau, cepgen::Particle::Status::FinalState);
  tau1.setMomentum(tau1_momentum);
  tau1.addMother(pho);
  ev.addParticle(tau1);
  cepgen::Particle tau2(
      cepgen::Particle::Role::CentralSystem, -(cepgen::spdgid_t)cepgen::PDG::tau, cepgen::Particle::Status::FinalState);
  tau2.setMomentum(tau2_momentum);
  tau2.addMother(pho);
  ev.addParticle(tau2);
  CG_LOG << ev;

  double weight = 1.;
  tauola->run(ev, weight);

  CG_TEST_SUMMARY;
}
