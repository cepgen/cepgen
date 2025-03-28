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

#include <cmath>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen::formfac {
  /// \cite Kelly:2004hm
  class Kelly final : public Parameterisation {
  public:
    explicit Kelly(const ParametersList& params)
        : Parameterisation(params),
          ae_(steer<std::vector<double> >("aE")),
          be_(steer<std::vector<double> >("bE")),
          am_(steer<std::vector<double> >("aM")),
          bm_(steer<std::vector<double> >("bM")) {}

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Kelly");
      desc.add("aE", std::vector{1., -0.24});
      desc.add("bE", std::vector{10.98, 12.82, 21.97});
      desc.add("aM", std::vector{1., 0.12});
      desc.add("bM", std::vector{10.97, 18.86, 6.55});
      return desc;
    }

  private:
    static double computeFF(double tau, const std::vector<double>& as, const std::vector<double>& bs) {
      double num{0.}, den{1.};
      for (size_t i = 0; i < as.size(); ++i)
        num += as.at(i) * std::pow(tau, i);
      for (size_t i = 0; i < bs.size(); ++i)
        den += bs.at(i) * std::pow(tau, i + 1);
      return num / den;
    }
    void eval() override {
      const auto ta = tau(q2_);
      setGEGM(computeFF(ta, ae_, be_), MU * computeFF(ta, am_, bm_));
    }

    const std::vector<double> ae_, be_, am_, bm_;
  };
}  // namespace cepgen::formfac
using cepgen::formfac::Kelly;
REGISTER_FORMFACTORS("Kelly", Kelly);
