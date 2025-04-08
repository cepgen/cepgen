/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2024  Laurent Forthomme
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

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();
  cepgen::initialise();

  auto tauola = cepgen::EventModifierFactory::get().build("tauola");
  if (!tauola) {
    CG_LOG << "Failed to retrieve the Tauola interface!";
    return -1;
  }
  tauola->initialise(cepgen::RunParameters());

  cepgen::Event ev;
  cepgen::Particle pho(cepgen::Particle::Role::CentralSystem, cepgen::PDG::photon, cepgen::Particle::Status::Resonance);
  ev.addParticle(pho);
  cepgen::Particle tau1(
      cepgen::Particle::Role::CentralSystem, +(cepgen::spdgid_t)cepgen::PDG::tau, cepgen::Particle::Status::FinalState);
  tau1.setMomentum(0., 0., 100.);
  tau1.addMother(pho);
  ev.addParticle(tau1);
  cepgen::Particle tau2(
      cepgen::Particle::Role::CentralSystem, -(cepgen::spdgid_t)cepgen::PDG::tau, cepgen::Particle::Status::FinalState);
  tau2.setMomentum(0., 0., -100.);
  tau2.addMother(pho);
  ev.addParticle(tau2);
  CG_LOG << ev;

  double weight = 1.;
  tauola->run(ev, weight);

  CG_TEST_SUMMARY;
}
