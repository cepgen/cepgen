/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6Interface.h"

namespace cepgen {
  namespace pythia6 {
    class AlphaEM final : public Coupling {
    public:
      explicit AlphaEM(const ParametersList& params) : Coupling(params) {
        mstu(101) = steer<int>("mode");
        paru(101) = steer<double>("fixedAlphaEM");
        paru(102) = steer<double>("sin2ThetaW");
        paru(103) = steer<double>("highQ2alphaEM");
        paru(104) = steer<double>("q2cut");
        paru(105) = steer<double>("gf");
      }

      inline static ParametersDescription description() {
        auto desc = cepgen::Coupling::description();
        desc.setDescription("Pythia6 modelling of alpha(EM) running");
        desc.add<int>("mode", mstu(101))
            .setDescription("procedure for alpha(EM) evaluation")
            .values()
            .allow(0, "fix at 'fixedAlphaEM'")
            .allow(1, "running accounting to fermion loops")
            .allow(2, "fix with low-high Q^2 splitting");
        desc.add<double>("fixedAlphaEM", paru(101))
            .setDescription("electromagnetic fine structure constant at vanishing mom.transfer");
        desc.add<double>("sin2ThetaW", paru(102)).setDescription("weak mixing angle of the standard electroweak model");
        desc.add<double>("highQ2alphaEM", paru(103))
            .setDescription("typical alpha(EM) in EW processes, intended for high-Q^2 for Z/W physics");
        desc.add<double>("q2cut", paru(104)).setDescription("dividing line between low- and high-Q^2 if mode=2");
        desc.add<double>("gf", paru(105)).setDescription("Fermi constant of weak interactions");
        return desc;
      }

      inline double operator()(double q) const override { return pyalem(q * q); }
    };
  }  // namespace pythia6
}  // namespace cepgen
using Pythia6AlphaEM = cepgen::pythia6::AlphaEM;
REGISTER_ALPHAEM_MODULE("pythia6", Pythia6AlphaEM);
