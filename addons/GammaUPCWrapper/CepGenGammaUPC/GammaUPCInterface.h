/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#ifndef CepGenAddOns_GammaUPCWrapper_GammaUPCInterface_h
#define CepGenAddOns_GammaUPCWrapper_GammaUPCInterface_h

#include <array>

namespace gammaUPC {
  struct to_collider_t {
    std::array<double, 2> ebeam, xbk, q2fact;
    std::array<int, 2> lpp;
  };
  struct ion_beam_t {
    int nuclearA_beam1, nuclearA_beam2, nuclearZ_beam1, nuclearZ_beam2;
  };
}  // namespace gammaUPC

extern "C" {
extern double __elasticphotonphotonflux_MOD_photonflux_proton_nob(double& x, double& gamma);
extern double __elasticphotonphotonflux_MOD_photonflux_nucleus_nob(double& x, double& gamma, double& z, double& ra);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pp(double& x1, double& x2, int& force_p_nohad1);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pA_hardsphere(double& x1, double& x2, int& force_p_nohad1);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pA_WoodsSaxon(double& x1, double& x2, int& force_p_nohad1);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_AB_hardsphere(double& x1, double& x2, int& force_p_nohad1);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_AB_WoodsSaxon(double& x1, double& x2, int& force_p_nohad1);
extern gammaUPC::to_collider_t to_collider_;
extern gammaUPC::ion_beam_t ion_beam_;
}

namespace gammaUPC {
  inline double photonFluxP(double x, double gamma) {
    return __elasticphotonphotonflux_MOD_photonflux_proton_nob(x, gamma);
  }

  inline double photonFluxA(double a, double gamma, double z, double ra) {
    return __elasticphotonphotonflux_MOD_photonflux_nucleus_nob(a, gamma, z, ra);
  }

  enum struct HeavyIonMode { HardSphere, WoodsSaxon };

  inline double twoPhotonFluxPP(double x1, double x2, int force_p_nohad1) {
    return __elasticphotonphotonflux_MOD_photonphotonflux_pp(x1, x2, force_p_nohad1);
  }
  inline double twoPhotonFluxPA(double x1, double x2, int force_p_nohad1, HeavyIonMode mode) {
    switch (mode) {
      case HeavyIonMode::HardSphere:
        return __elasticphotonphotonflux_MOD_photonphotonflux_pA_hardsphere(x1, x2, force_p_nohad1);
      case HeavyIonMode::WoodsSaxon:
        return __elasticphotonphotonflux_MOD_photonphotonflux_pA_WoodsSaxon(x1, x2, force_p_nohad1);
      default:
        return -1.;
    }
  }
  inline double twoPhotonFluxAB(double x1, double x2, int force_p_nohad1, HeavyIonMode mode) {
    switch (mode) {
      case HeavyIonMode::HardSphere:
        return __elasticphotonphotonflux_MOD_photonphotonflux_AB_hardsphere(x1, x2, force_p_nohad1);
      case HeavyIonMode::WoodsSaxon:
        return __elasticphotonphotonflux_MOD_photonphotonflux_AB_WoodsSaxon(x1, x2, force_p_nohad1);
      default:
        return -1.;
    }
  }
}  // namespace gammaUPC

#endif