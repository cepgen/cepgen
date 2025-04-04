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

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

namespace cepgen::formfac {
  class PointLike : public Parameterisation {
  public:
    explicit PointLike(const ParametersList& params, double fe, double fm)
        : Parameterisation(params), fe_(fe), fm_(fm) {}

  private:
    void eval() override { setFEFM(fe_, fm_); }
    const double fe_, fm_;
  };

  struct PointLikeScalar final : PointLike {
    explicit PointLikeScalar(const ParametersList& params) : PointLike(params, 1., 0.) {}
    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Point-like scalar");
      return desc;
    }
  };

  struct PointLikeFermion final : PointLike {
    explicit PointLikeFermion(const ParametersList& params) : PointLike(params, 1., 1.) {}
    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Point-like fermion");
      return desc;
    }
  };
}  // namespace cepgen::formfac
using cepgen::formfac::PointLikeFermion;
using cepgen::formfac::PointLikeScalar;
REGISTER_FORMFACTORS("PointLikeScalar", PointLikeScalar);
REGISTER_FORMFACTORS("PointLikeFermion", PointLikeFermion);
