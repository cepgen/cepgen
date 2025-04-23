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

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGenMadGraph/Interface.h"
#include "CepGenMadGraph/Process.h"
#include "CepGenMadGraph/ProcessBuilder.h"

using namespace cepgen;
using namespace std::string_literals;

namespace cepgen::mg5amc {
  class GeneralProcessBuilder final : public ProcessBuilder {
  public:
    explicit GeneralProcessBuilder(const ParametersList& params, bool load_library = true)
        : ProcessBuilder(params, load_library) {
      setCentral(mg5_proc_->centralSystem());
    }

    proc::ProcessPtr clone() const override { return std::make_unique<GeneralProcessBuilder>(parameters(), false); }

    void prepareFactorisedPhaseSpace() override {
      if (const auto psgen_partons = phase_space_generator_->partons();
          mg5_proc_->intermediatePartons() != psgen_partons)
        throw CG_FATAL("mg5amc:GeneralProcessBuilder")
            << "MadGraph unpacked process incoming state (" << mg5_proc_->intermediatePartons() << ") "
            << "is incompatible with user-steered incoming fluxes particles (" << psgen_partons << ").";
      if (const auto params_card = steer<std::string>("parametersCard"); !params_card.empty()) {
        CG_INFO("mg5amc:GeneralProcessBuilder")
            << "Preparing process kinematics for card at \"" << params_card << "\".";
        const auto unsteered_pcard = Interface::extractParamCardParameters(utils::readFile(params_card));
        CG_DEBUG("mg5amc:GeneralProcessBuilder") << "Unsteered parameters card:\n" << unsteered_pcard;
        if (const auto mod_params = steer<ParametersList>("modelParameters"); !mod_params.empty()) {
          const auto steered_pcard = unsteered_pcard.steer(mod_params);
          CG_DEBUG("mg5amc:GeneralProcessBuilder") << "User-steered parameters:" << mod_params << "\n"
                                                   << "Steered parameters card:\n"
                                                   << steered_pcard;
          std::ofstream params_card_steered(params_card);
          params_card_steered << Interface::generateParamCard(steered_pcard);
          params_card_steered.close();
        }
        mg5_proc_->initialise(params_card);
      }
    }
    double computeFactorisedMatrixElement() override {
      if (!mg5_proc_)
        CG_FATAL("mg5amc:GeneralProcessBuilder:eval") << "Process not properly linked!";
      if (!kinematics().cuts().initial.contain(event()(Particle::Role::Parton1)) ||
          !kinematics().cuts().initial.contain(event()(Particle::Role::Parton2)))
        return 0.;
      if (!kinematics().cuts().central.contain(event()(Particle::Role::CentralSystem)))
        return 0.;

      CG_DEBUG_LOOP("mg5amc:GeneralProcessBuilder:eval")
          << "Particles content:\n"
          << "incoming: " << q1() << " (m=" << q1().mass() << "), " << q2() << " (m=" << q2().mass() << ")\n"
          << "outgoing: " << pc(0) << " (m=" << pc(0).mass() << "), " << pc(1) << " (m=" << pc(1).mass() << ").";
      mg5_proc_->setMomentum(0, q1());   // first incoming parton
      mg5_proc_->setMomentum(1, q2());   // second incoming parton
      mg5_proc_->setMomentum(2, pc(0));  // first outgoing central particle
      mg5_proc_->setMomentum(3, pc(1));  // second outgoing central particle
      if (const auto weight = mg5_proc_->eval(); utils::positive(weight))
        return weight * std::pow(shat(), -2);
      return 0.;
    }
  };
}  // namespace cepgen::mg5amc
using MadGraphGeneralProcessBuilder = mg5amc::GeneralProcessBuilder;
REGISTER_PROCESS("mg5_aMC", MadGraphGeneralProcessBuilder);
