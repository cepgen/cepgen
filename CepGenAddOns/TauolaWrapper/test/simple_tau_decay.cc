/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include "CepGen/Core/EventModifier.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/EventModifierFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"
#include "CepGenAddOns/TauolaWrapper/PhotosTauolaInterface.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  auto tauola = cepgen::EventModifierFactory::get().build("tauola");
  if (!tauola) {
    CG_LOG << "Failed to retrieve the Tauola interface!";
    return -1;
  }
  tauola->init();

  cepgen::Event ev;
  cepgen::Particle tau(cepgen::Particle::Role::CentralSystem, cepgen::PDG::tau, cepgen::Particle::Status::FinalState);
  ev.addParticle(tau);

  double weight = 1.;
  tauola->run(ev, weight);

  return 0;
}
