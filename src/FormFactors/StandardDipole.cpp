/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/PDG.h"

namespace cepgen::formfac {
  class StandardDipole : public Parameterisation {
  public:
    explicit StandardDipole(const ParametersList& params)
        : Parameterisation(params), inv_sq_scale_param_(1. / steer<double>("scale")) {}

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("Standard dipole");
      desc.add<pdgid_t>("pdgId", PDG::proton);
      desc.add("scale", 0.71).setDescription("scaling (in GeV^2) (0.71 for r_p = 0.81 fm, 0.66 for r_p = 0.84 fm)");
      return desc;
    }

  protected:
    void eval() override {
      const auto ge = pow(1. + q2_ * inv_sq_scale_param_, -2.);
      setGEGM(ge, MU * ge);
    }

  private:
    const double inv_sq_scale_param_;
  };

  class HeavyIonDipole final : public StandardDipole {
  public:
    explicit HeavyIonDipole(const ParametersList& params)
        : StandardDipole(params),
          hi_(HeavyIon::fromPdgId(pdg_id_)),
          a2_(std::pow(hi_.radius() / constants::GEVM1_TO_M, 2)),
          a02_(std::pow(HeavyIon::proton().radius() / constants::GEVM1_TO_M, 2)) {}

    static ParametersDescription description() {
      auto desc = StandardDipole::description();
      desc.setDescription("Heavy ion dipole");
      desc.addAs<pdgid_t, HeavyIon>("pdgId", HeavyIon::Pb());
      return desc;
    }

  private:
    void eval() override {
      if (hi_ == HeavyIon::proton()) {
        StandardDipole::eval();
        return;
      }
      const auto qr2 = q2_ * a2_;
      if (static_cast<short>(hi_.Z) <= static_cast<short>(Element::C)) {  // Gaussian form factor for light nuclei
        const auto ge = std::exp(-qr2 / 6.);
        setGEGM(ge, MU * ge);
        return;
      }
      const auto qr = std::sqrt(qr2), inv_qr = 1. / qr;
      const auto sph = (std::sin(qr) - qr * std::cos(qr)) * 3. * inv_qr * inv_qr;
      const auto ge = sph / (1. + q2_ * a02_);
      setGEGM(ge, MU * ge);
    }
    const HeavyIon hi_;
    const double a2_, a02_;
  };
}  // namespace cepgen::formfac
using cepgen::formfac::HeavyIonDipole;
using cepgen::formfac::StandardDipole;
REGISTER_FORMFACTORS(cepgen::formfac::gFFStandardDipoleHandler, StandardDipole);
REGISTER_FORMFACTORS("HeavyIonDipole", HeavyIonDipole);
