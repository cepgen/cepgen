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

#include "CepGenGammaUPC/GammaUPCInterface.h"

#include <CepGen/Core/Exception.h>

extern "C" {
extern double __elasticphotonphotonflux_MOD_photonflux_proton_nob(double& x, double& gamma);
extern double __elasticphotonphotonflux_MOD_photonflux_nucleus_nob(double& x, double& gamma, double& z, double& ra);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pp(double& x1, double& x2, int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pA_hardsphere(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pA_WoodsSaxon(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_AB_hardsphere(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_AB_WoodsSaxon(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
}

namespace gammaUPC {
  double photonFluxP(double x, double gamma) {
    const auto flux = __elasticphotonphotonflux_MOD_photonflux_proton_nob(x, gamma);
    CG_DEBUG_LOOP("gammaUPC:photonFluxP")
        << "Single-photon flux from proton (gamma=" << gamma << ", no b-effect): " << flux << ".";
    return flux;
  }

  double photonFluxA(double a, double gamma, double z, double ra) {
    const auto flux = __elasticphotonphotonflux_MOD_photonflux_nucleus_nob(a, gamma, z, ra);
    CG_DEBUG_LOOP("gammaUPC:photonFluxA") << "Single-photon flux from HI{A=" << a << ", Z=" << z << ", gamma=" << gamma
                                          << ", rA=" << ra << "} (no b-effect): " << flux << ".";
    return flux;
  }

  double twoPhotonFluxPP(double x1, double x2, int force_p_no_hadronisation) {
    const auto flux = __elasticphotonphotonflux_MOD_photonphotonflux_pp(x1, x2, force_p_no_hadronisation);
    CG_DEBUG_LOOP("gammaUPC:twoPhotonFluxPP")
        << "Two-photon flux from proton-proton for x1=" << x1 << ", x2=" << x2
        << " [non-had.: " << static_cast<bool>(force_p_no_hadronisation) << "]: " << flux << ".";
    return flux;
  }

  double twoPhotonFluxPA(double x1, double x2, int force_p_no_hadronisation, HeavyIonMode mode) {
    switch (mode) {
      case HeavyIonMode::HardSphere: {
        const auto flux =
            __elasticphotonphotonflux_MOD_photonphotonflux_pA_hardsphere(x1, x2, force_p_no_hadronisation);
        CG_DEBUG_LOOP("gammaUPC:twoPhotonFluxPA")
            << "Two-photon flux from proton-HI with hard sphere modelling for x1=" << x1 << ", x2=" << x2
            << " [non-had.: " << static_cast<bool>(force_p_no_hadronisation) << "]: " << flux << ".";
        return flux;
      }
      case HeavyIonMode::WoodsSaxon: {
        const auto flux =
            __elasticphotonphotonflux_MOD_photonphotonflux_pA_WoodsSaxon(x1, x2, force_p_no_hadronisation);
        CG_DEBUG_LOOP("gammaUPC:twoPhotonFluxPA")
            << "Two-photon flux from proton-HI with Woods-Saxon modelling for x1=" << x1 << ", x2=" << x2
            << " [non-had.: " << static_cast<bool>(force_p_no_hadronisation) << "]: " << flux << ".";
        return flux;
      }
      default:
        throw CG_FATAL("gammaUPC:twoPhotonFluxPA")
            << "Invalid mode for the heavy ion modelling: " << static_cast<int>(mode) << ".";
    }
  }

  double twoPhotonFluxAB(double x1, double x2, int force_p_no_hadronisation, HeavyIonMode mode) {
    switch (mode) {
      case HeavyIonMode::HardSphere: {
        const auto flux =
            __elasticphotonphotonflux_MOD_photonphotonflux_AB_hardsphere(x1, x2, force_p_no_hadronisation);
        CG_DEBUG_LOOP("gammaUPC:twoPhotonFluxAB")
            << "Two-photon flux from HI-HI with hard sphere modelling for x1=" << x1 << ", x2=" << x2
            << " [non-had.: " << static_cast<bool>(force_p_no_hadronisation) << "]: " << flux << ".";
        return flux;
      }
      case HeavyIonMode::WoodsSaxon: {
        const auto flux =
            __elasticphotonphotonflux_MOD_photonphotonflux_AB_WoodsSaxon(x1, x2, force_p_no_hadronisation);
        CG_DEBUG_LOOP("gammaUPC:twoPhotonFluxAB")
            << "Two-photon flux from HI-HI with Woods-Saxon modelling for x1=" << x1 << ", x2=" << x2
            << " [non-had.: " << static_cast<bool>(force_p_no_hadronisation) << "]: " << flux << ".";
        return flux;
      }
      default:
        throw CG_FATAL("gammaUPC:twoPhotonFluxAB")
            << "Invalid mode for the heavy ion modelling: " << static_cast<int>(mode) << ".";
    }
  }
}  // namespace gammaUPC
