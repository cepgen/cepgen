/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

namespace cepgen {
  /// A perturbative PDF-oriented \f$\alpha_S(Q^2)\f$ evaluator
  class AlphaSLHAPDFAnalytic final : public Coupling {
  public:
    explicit AlphaSLHAPDFAnalytic(const ParametersList& params) : Coupling(params), ana_(new LHAPDF::AlphaS_Analytic) {
      ana_->setOrderQCD(steer<int>("order"));
      for (int i = 1; i <= 6; ++i)  // set all quarks masses for evolution
        ana_->setQuarkMass(i, PDG::get().mass(i));
      // set gradients for evolution
      size_t i = 3;
      for (const auto& lambda : steer<std::vector<double> >("lambdas"))
        ana_->setLambda(i++, lambda);
    }

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("Analytic LHAPDF perturb.algo.");
      desc.add<std::string>("pdfSet", "cteq66");
      desc.add<int>("order", 4).setDescription("QCD order");
      desc.add<std::vector<double> >("lambdas", {0.339, 0.296, 0.213});
      return desc;
    }

    double operator()(double q) const override { return ana_->alphasQ(q); }

  private:
    std::unique_ptr<LHAPDF::AlphaS_Analytic> ana_;
  };
}  // namespace cepgen

REGISTER_ALPHAS_MODULE("lhapdfAnalytic", AlphaSLHAPDFAnalytic);

#endif
