/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2024  Laurent Forthomme
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

#include <APFEL/APFEL.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen {
  namespace apfel {
    class AlphaS final : public cepgen::Coupling {
    public:
      explicit AlphaS(const ParametersList& params) : cepgen::Coupling(params), q_range_(steer<Limits>("qrange")) {
        APFEL::SetPerturbativeOrder(steer<int>("order"));
        APFEL::SetPoleMasses(PDG::get().mass(4), PDG::get().mass(5), PDG::get().mass(6));
        APFEL::InitializeAPFEL();
        APFEL::EvolveAPFEL(q_range_.min(), q_range_.max());
        if (steer<bool>("checkAPFEL") && !APFEL::CheckAPFEL())
          throw CG_FATAL("apfel:AlphaS") << "Something is wrong with your APFEL configuration.";
      }
      static ParametersDescription description() {
        auto desc = cepgen::Coupling::description();
        desc.setDescription("APFEL alpha(S) evolution algorithm");
        desc.add<bool>("checkAPFEL", false).setDescription("perform full check of APFEL configuration");
        desc.add<int>("order", 2).setDescription("QCD perturbative evolution order");
        desc.add<Limits>("qrange", {1., 1.e4}).setDescription("Q range reachable for evolution (in GeV)");
        return desc;
      }

      double operator()(double q) const override {
        if (!q_range_.contains(q))
          CG_WARNING("apfel:AlphaS:get") << "q = " << q << " outside the evolution range" << q_range_ << ".";
        return APFEL::AlphaQCD(q);
      }

    private:
      const Limits q_range_;
    };
  }  // namespace apfel
}  // namespace cepgen
using AlphaSAPFEL = cepgen::apfel::AlphaS;
REGISTER_ALPHAS_MODULE("apfel", AlphaSAPFEL);
