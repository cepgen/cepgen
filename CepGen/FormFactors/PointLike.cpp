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

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen {
  namespace formfac {
    class PointLike : public Parameterisation {
    public:
      explicit PointLike(const ParametersList& params, const FormFactors& ff)
          : Parameterisation(params), trivial_(ff) {}

    private:
      FormFactors compute(double /*q2*/) override { return trivial_; }
      const FormFactors trivial_;
    };

    struct PointLikeScalar final : public PointLike {
      explicit PointLikeScalar(const ParametersList& params) : PointLike(params, FormFactors{1., 0., 0., 0.}) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Point-like scalar");
        return desc;
      }
    };

    struct PointLikeFermion final : public PointLike {
      explicit PointLikeFermion(const ParametersList& params) : PointLike(params, FormFactors{1., 1., 0., 0.}) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Point-like fermion");
        return desc;
      }
    };
  }  // namespace formfac
}  // namespace cepgen

REGISTER_FORMFACTORS("PointLikeScalar", PointLikeScalar)
REGISTER_FORMFACTORS("PointLikeFermion", PointLikeFermion)
