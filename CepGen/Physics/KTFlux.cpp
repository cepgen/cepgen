/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  const double KTFluxParameters::kMinKTFlux = 1.e-20;

  double ktFlux(const KTFlux& type,
                double x,
                double kt2,
                formfac::Parameterisation& ff,
                strfun::Parameterisation& sf,
                double mi2,
                double mf2) {
    switch (type) {
      case KTFlux::P_Photon_Elastic:
      case KTFlux::P_Photon_Elastic_Budnev: {
        const double x2 = x * x;
        const double q2min = x2 * mi2 / (1. - x), q2 = q2min + kt2 / (1. - x);
        const double qnorm = 1. - q2min / q2;
        const auto& formfac = ff(mode::Beam::ProtonElastic, q2);
        if (type == KTFlux::P_Photon_Elastic) {
          const double f_aux = formfac.FE * qnorm * qnorm;
          return constants::ALPHA_EM * M_1_PI / q2 * f_aux;
        } else {
          const double f_D = formfac.FE * (1. - x) * qnorm;
          const double f_C = formfac.FM;
          return constants::ALPHA_EM * M_1_PI * (1. - x) / q2 * (f_D + 0.5 * x2 * f_C);
        }
      } break;
      case KTFlux::P_Photon_Inelastic:
      case KTFlux::P_Photon_Inelastic_Budnev: {
        const double x2 = x * x;
        const double q2min = (x * (mf2 - mi2) + x2 * mi2) / (1. - x);
        const double q2 = q2min + kt2 / (1. - x);
        const double qnorm = 1. - q2min / q2;
        //--- proton structure functions
        const double denom = 1. / (q2 + mf2 - mi2);
        const double xbj = denom * q2;
        //--- proton structure functions
        if (type == KTFlux::P_Photon_Inelastic) {
          const double f_aux = sf.F2(xbj, q2) * denom * qnorm * qnorm;
          return constants::ALPHA_EM * M_1_PI * (1. - x) / q2 * f_aux;
        } else {
          const double f_D = sf.F2(xbj, q2) * denom * (1. - x) * qnorm;
          const double f_C = sf.F1(xbj, q2) * 2. / q2;
          return constants::ALPHA_EM * M_1_PI * (1. - x) / q2 * (f_D + 0.5 * x2 * f_C);
        }
      } break;
      case KTFlux::P_Gluon_KMR: {
        return kmr::GluonGrid::get()(x, kt2, mf2);
      } break;
      default:
        throw CG_FATAL("KTFlux") << "Invalid flux type: " << type;
    }
  }

  double ktFlux(const KTFlux& type, double x, double kt2, const HeavyIon& hi) {
    const double& mp = PDG::get().mass(PDG::proton);
    double flux = 0.;
    switch (type) {
      case KTFlux::HI_Photon_Elastic: {
        const double r_a = 1.1 * cbrt(hi.A), a0 = 0.7, m_a = hi.A * mp;
        const double q2_ela = (kt2 + x * x * m_a * m_a) / (1. - x), cons = sqrt(q2_ela) / 0.1973;
        const double tau = cons * r_a, tau1 = cons * a0;
        // "Realistic nuclear form-factor" as used in STARLIGHT
        const double ff1 = 3. * (sin(tau) - tau * cos(tau)) / pow(tau + 1.e-10, 3);
        const double ff2 = 1. / (1. + tau1 * tau1);
        const double ela1 = pow(kt2 / (kt2 + x * x * m_a * m_a), 2);
        const double ela2 = pow(ff1 * ff2, 2) /*, ela3 = 1.-( q2_ela-kt2 )/q2_ela*/;
        const unsigned int z = (unsigned short)hi.Z;
        flux = constants::ALPHA_EM * M_1_PI * z * z * ela1 * ela2 / q2_ela;
      } break;
      default:
        throw CG_FATAL("KTFlux") << "Invalid flux type: " << type;
    }
    if (flux < KTFluxParameters::kMinKTFlux)
      return 0.;
    return flux;
  }

  std::ostream& operator<<(std::ostream& os, const KTFlux& type) {
    switch (type) {
      case KTFlux::P_Photon_Elastic:
        return os << "elastic photon from proton";
      case KTFlux::P_Photon_Elastic_Budnev:
        return os << "elastic photon from proton (Budnev)";
      case KTFlux::P_Photon_Inelastic:
        return os << "inelastic photon from proton";
      case KTFlux::P_Photon_Inelastic_Budnev:
        return os << "inelastic photon from proton (Budnev)";
      case KTFlux::P_Gluon_KMR:
        return os << "elastic gluon from proton (KMR)";
      case KTFlux::HI_Photon_Elastic:
        return os << "elastic photon from HI";
      case KTFlux::invalid:
      default:
        return os << "unrecognised flux (" << (int)type << ")";
    }
  }
}  // namespace cepgen
