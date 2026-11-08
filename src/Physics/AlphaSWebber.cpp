/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Utils/Message.h"

using namespace cepgen;

/// Simple parameterisation of the QCD running coupling at low scales
/// \cite Webber:1998um
class AlphaSWebber final : public Coupling {
public:
  explicit AlphaSWebber(const ParametersList& params)
      : Coupling(params),
        num_colours_(steer<int>("Nc")),
        num_flavours_(steer<int>("nf")),
        lambda_(steer<double>("Lambda")),
        beta0_((11. * num_colours_ - 2. * num_flavours_) / 3.),
        prefactor_(4. * M_PI / beta0_) {
    CG_INFO("AlphaSWebber:init") << "Webber et al. alpha(S) evolution algorithm initialised with parameters:\n\t"
                                 << "Nc: " << num_colours_ << ", nf: " << num_flavours_ << " -> beta0: " << beta0_
                                 << ", Lambda: " << lambda_ << ".";
  }

  static ParametersDescription description() {
    auto desc = Coupling::description();
    desc.setDescription("Webber alpha(S) evolution algorithm");
    desc.add("Nc", 3).setDescription("number of colours considered");
    desc.add("nf", 3).setDescription("number of fermion flavours considered");
    desc.add("Lambda", 0.25).setDescription("evolution scale (in GeV)");
    return desc;
  }

  double operator()(double q) const override {
    const auto mun = q * q / lambda_ / lambda_;
    return prefactor_ * (1. / std::log(mun) + 125. * (1. + 4. * mun) / (1. - mun) / std::pow(4. + mun, 4));
  }

private:
  const int num_colours_, num_flavours_;
  const double lambda_, beta0_;
  const double prefactor_;
};
REGISTER_ALPHAS_MODULE("webber", AlphaSWebber);
