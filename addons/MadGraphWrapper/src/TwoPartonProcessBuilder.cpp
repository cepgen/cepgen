/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGenMadGraph/Process.h"
#include "CepGenMadGraph/ProcessBuilder.h"

namespace cepgen::mg5amc {
  class TwoPartonProcessBuilder final : public ProcessBuilder {
  public:
    explicit TwoPartonProcessBuilder(const ParametersList& params, bool load_library = true)
        : ProcessBuilder(params, load_library) {
      setCentral(process().centralSystem());
    }

    proc::ProcessPtr clone() const override { return std::make_unique<TwoPartonProcessBuilder>(parameters(), false); }

    void prepareFactorisedPhaseSpace() override {
      if (const auto psgen_partons = phase_space_generator_->partons();
          process().intermediatePartons() != psgen_partons)
        throw CG_FATAL("mg5amc:TwoPartonProcessBuilder")
            << "MadGraph unpacked process incoming state (" << process().intermediatePartons() << ") "
            << "is incompatible with user-steered incoming fluxes particles (" << psgen_partons << ").";
      prepareSteeringCard();
    }
    double computeFactorisedMatrixElement() override {
      if (!kinematics().cuts().initial.contain(event()(Particle::Role::Parton1)) ||
          !kinematics().cuts().initial.contain(event()(Particle::Role::Parton2)))
        return 0.;
      if (!kinematics().cuts().central.contain(event()(Particle::Role::CentralSystem)))
        return 0.;

      CG_DEBUG_LOOP("mg5amc:TwoPartonProcessBuilder:eval")
          << "Particles content:\n"
          << "incoming: " << q1() << " (m=" << q1().mass() << "), " << q2() << " (m=" << q2().mass() << ")\n"
          << "outgoing: " << pc(0) << " (m=" << pc(0).mass() << "), " << pc(1) << " (m=" << pc(1).mass() << ").";
      process().setMomentum(0, q1());  // first incoming parton
      process().setMomentum(1, q2());  // second incoming parton
      for (size_t i = 0; i < event()(Particle::Role::CentralSystem).size(); ++i)
        process().setMomentum(2 + i, pc(i));  // outgoing central particles
      if (const auto weight = process().eval(); utils::positive(weight))
        return weight * std::pow(shat(), -2);
      return 0.;
    }
  };
}  // namespace cepgen::mg5amc
using MadGraphTwoPartonProcessBuilder = cepgen::mg5amc::TwoPartonProcessBuilder;
REGISTER_PROCESS("mg5_aMC", MadGraphTwoPartonProcessBuilder);
