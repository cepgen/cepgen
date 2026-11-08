/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#include <LHAPDF/LHAPDF.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"

#if defined LHAPDF_MAJOR_VERSION && LHAPDF_MAJOR_VERSION == 6

namespace cepgen::lhapdf {
  /// A perturbative PDF-oriented \f$\alpha_S(Q^2)\f$ evaluator
  class AlphaSODE final : public Coupling {
  public:
    explicit AlphaSODE(const ParametersList& params) : Coupling(params), ode_(new LHAPDF::AlphaS_ODE) {
      ode_->setOrderQCD(steer<int>("order"));
      ode_->setAlphaSMZ(steer<double>("alphaSMZ"));
      ode_->setMZ(PDG::get().mass(23));  // set Z mass
      for (int i = 1; i <= 6; ++i)       // set all quarks masses
        ode_->setQuarkMass(i, PDG::get().mass(i));
    }

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("ODE LHAPDF evol.algo.");
      desc.add("order", 5).setDescription("QCD order");
      desc.add("alphaSMZ", 0.118);
      return desc;
    }

    double operator()(double q) const override { return ode_->alphasQ(q); }

  private:
    const std::unique_ptr<LHAPDF::AlphaS_ODE> ode_;
  };
}  // namespace cepgen::lhapdf
using AlphaS_LHAPDF_ODE = cepgen::lhapdf::AlphaSODE;
REGISTER_ALPHAS_MODULE("lhapdfODE", AlphaS_LHAPDF_ODE);

#endif
