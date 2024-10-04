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
#include "CepGenPythia6/Pythia6Interface.h"

namespace cepgen::pythia6 {
  class AlphaS final : public Coupling {
  public:
    explicit AlphaS(const ParametersList& params) : Coupling(params) {
      mstu(111) = steer<int>("order");
      mstu(112) = steer<int>("nf");
      mstu(113) = steer<int>("minNf");
      mstu(114) = steer<int>("maxNf");
      mstu(115) = steer<int>("singularityTreatment");
      paru(111) = steer<double>("fixedAlphaS");
      paru(112) = steer<double>("Lambda");
      paru(113) = steer<double>("flavourThreshold");
      paru(114) = steer<double>("minQ2");
      paru(115) = steer<double>("maxQ2");
    }

    inline static ParametersDescription description() {
      auto desc = cepgen::Coupling::description();
      desc.setDescription("Pythia6 modelling of alpha(S) running");
      desc.add<int>("order", mstu(111))
          .setDescription("order of alpha(S) evaluation")
          .allow(0, "fixed at 'fixedAlphaS'")
          .allow(1, "1st order running")
          .allow(2, "2nd order running");
      desc.add<int>("nf", mstu(112)).setDescription("nominal number of ﬂavours assumed in alpha(s) expression");
      desc.add<int>("minNf", mstu(113))
          .setDescription("minimum number of ﬂavours that may be assumed in alpha(S) expression");
      desc.add<int>("maxNf", mstu(114))
          .setDescription("minimum number of ﬂavours that may be assumed in alpha(S) expression");
      desc.add<int>("singularityTreatment", mstu(115))
          .setDescription("treatment of alpha(S) singularities for Q^2->0")
          .allow(0, "allow divergence")
          .allow(1, "log-softening")
          .allow(2, "freeze under Q^2 transition value");
      desc.add<double>("fixedAlphaS", paru(111))
          .setDescription(
              "fix alpha(S) value assumed when order=0 (and also in parton showers when alpha(S) is assumed fix "
              "there)");
      desc.add<double>("Lambda", paru(112)).setDescription("Lambda value used in running");
      desc.add<double>("flavourThreshold", paru(113))
          .setDescription(
              "flavour threshold, for the effective number of flavours 'nf' to use (='flavourThreshold'*m_q^2)");
      desc.add<double>("minQ2", paru(114)).setDescription("Q^2 value below which alpha(S) is assumed constant");
      desc.add<double>("maxQ2", paru(115)).setDescription("maximum alpha(S) value computable");
      return desc;
    }

    inline double operator()(double q) const override { return pyalps(q * q); }
  };
}  // namespace cepgen::pythia6
using Pythia6AlphaS = cepgen::pythia6::AlphaS;
REGISTER_ALPHAS_MODULE("pythia6", Pythia6AlphaS);
