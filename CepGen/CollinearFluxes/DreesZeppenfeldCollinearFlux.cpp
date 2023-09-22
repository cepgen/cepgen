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

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  /// Virtuality-dependent Drees-Zeppenfeld photon flux
  /// \note Corresponds to  PDF:Proton2gammaSet=2 in Pythia 8
  /// \cite Drees:1988pp
  class DreesZeppenfeldCollinearFlux : public CollinearFlux {
  public:
    explicit DreesZeppenfeldCollinearFlux(const ParametersList& params)
        : CollinearFlux(params), scale_(steer<double>("scale")) {}

    static ParametersDescription description() {
      auto desc = CollinearFlux::description();
      desc.setDescription("Drees-Zeppenfeld Q^2-dependent flux");
      desc.add<double>("scale", 0.71);
      return desc;
    }

    bool fragmenting() const override final { return true; }
    pdgid_t partonPdgId() const override final { return PDG::photon; }

    double fluxQ2(double x, double q2) const override {
      if (!x_range_.contains(x, true))
        return 0.;
      const auto fq4 = std::pow(1 + q2 / scale_, -4);  // Q^2-dependent form factor
      return prefactor_ * 0.5 * (1. + std::pow(1. - x, 2)) / q2 * fq4;
    }

  protected:
    double mass2() const override { return mp2_; }
    const double scale_;
  };
}  // namespace cepgen
REGISTER_COLLINEAR_FLUX("DreesZeppenfeld", DreesZeppenfeldCollinearFlux);
