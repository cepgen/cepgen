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

using namespace std::string_literals;

namespace cepgen::lhapdf {
  /// A perturbative PDF-oriented \f$\alpha_S(Q^2)\f$ evaluator
  class AlphaSAnalytic final : public Coupling {
  public:
    explicit AlphaSAnalytic(const ParametersList& params)
        : Coupling(params), alphas_analytic_(new LHAPDF::AlphaS_Analytic) {
      alphas_analytic_->setOrderQCD(steer<int>("order"));
      for (int i = 1; i <= 6; ++i)  // set all quarks masses for evolution
        alphas_analytic_->setQuarkMass(i, PDG::get().mass(i));
      // set gradients for evolution
      size_t i = 3;
      for (const auto& lambda : steer<std::vector<double> >("lambdas"))
        alphas_analytic_->setLambda(i++, lambda);
    }

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("Analytic LHAPDF perturb.algo.");
      desc.add("pdfSet", "cteq66"s);
      desc.add("order", 4).setDescription("QCD order");
      desc.add("lambdas", std::vector{0.339, 0.296, 0.213});
      return desc;
    }

    double operator()(double q) const override { return alphas_analytic_->alphasQ(q); }

  private:
    const std::unique_ptr<LHAPDF::AlphaS_Analytic> alphas_analytic_;
  };
}  // namespace cepgen::lhapdf
using AlphaS_LHAPDFAnalytic = cepgen::lhapdf::AlphaSAnalytic;
REGISTER_ALPHAS_MODULE("lhapdfAnalytic", AlphaS_LHAPDFAnalytic);

#endif
