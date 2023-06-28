/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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
#include "CepGen/Utils/Message.h"

namespace cepgen {
  namespace formfac {
    class StandardDipole : public Parameterisation {
    public:
      explicit StandardDipole(const ParametersList& params)
          : Parameterisation(params), inv_sq_scale_param_(1. / steer<double>("scale")) {}

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("Standard dipole");
        desc.add<pdgid_t>("pdgId", PDG::proton);
        desc.add<double>("scale", 0.71)
            .setDescription("scaling (in GeV^2) (0.71 for r_p = 0.81 fm, 0.66 for r_p = 0.84 fm)");
        return desc;
      }

    protected:
      FormFactors compute(double q2) override {
        FormFactors out;
        out.GE = pow(1. + q2 * inv_sq_scale_param_, -2.);
        out.GM = MU * out.GE;
        return out;
      }

    private:
      const double inv_sq_scale_param_;
    };

    class HeavyIonDipole final : public StandardDipole {
    public:
      explicit HeavyIonDipole(const ParametersList& params)
          : StandardDipole(params),
            hi_(HeavyIon::fromPdgId(pdg_id_)),
            a_(hi_.radius() / (constants::GEVM1_TO_M * 1e15)),
            a0_(HeavyIon::proton().radius() / (constants::GEVM1_TO_M * 1e15)),  // [fm -> GeV]
            a02_(a0_ * a0_) {}

      static ParametersDescription description() {
        auto desc = StandardDipole::description();
        desc.setDescription("Heavy ion dipole");
        desc.addAs<pdgid_t, HeavyIon>("pdgId", HeavyIon::Pb());
        return desc;
      }

    private:
      FormFactors compute(double q2) override {
        if (hi_ == HeavyIon::proton())
          return StandardDipole::compute(q2);
        if ((short)hi_.Z < 7) {  // Gaussian form factor for light nuclei
          FormFactors out;
          out.GE = exp(-a_ * a_ * q2 / 6.);
          out.GM = MU * out.GE;
          return out;
        }
        const double qr = sqrt(q2) * a_, inv_qr = 1. / qr;
        const double sph = (sin(qr) - qr * cos(qr)) * 3. * inv_qr * inv_qr * inv_qr;
        FormFactors out;
        out.GE = sph / (1. + q2 * a02_);
        out.GM = MU * out.GE;
        return out;
      }
      const HeavyIon hi_;
      const double a_, a0_, a02_;
    };
  }  // namespace formfac
}  // namespace cepgen
typedef cepgen::formfac::StandardDipole DipoleFF;
typedef cepgen::formfac::HeavyIonDipole HIDipoleFF;
REGISTER_FORMFACTORS(cepgen::formfac::gFFStandardDipoleHandler, DipoleFF);
REGISTER_FORMFACTORS("HeavyIonDipole", HIDipoleFF);
