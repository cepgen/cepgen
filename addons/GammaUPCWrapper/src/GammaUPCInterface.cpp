/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2026  Laurent Forthomme
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

#include <CepGen/Core/Exception.h>

#include <ostream>

#include "CepGenGammaUPC/GammaUPCInterface.h"

// single-photon (no impact parameter dependence) fluxes: names depend on the gammaUPC version -> resolve when loading
// (v1.0 standalone: photonflux_*_nob; v1.7, as shipped with MadGraph5_aMC@NLO: photonflux_*_edff_nob)
extern "C" {
extern double __elasticphotonphotonflux_MOD_photonflux_proton_nob(double& x, double& gamma) __attribute__((weak));
extern double __elasticphotonphotonflux_MOD_photonflux_proton_edff_nob(double& x, double& gamma) __attribute__((weak));
extern double __elasticphotonphotonflux_MOD_photonflux_proton_chff_nob(double& x, double& gamma, int* correction)
    __attribute__((weak));  // v1.7, optional LOGICAL argument (nullptr if absent)
extern double __elasticphotonphotonflux_MOD_photonflux_nucleus_nob(double& x, double& gamma, double& z, double& ra)
    __attribute__((weak));
extern double __elasticphotonphotonflux_MOD_photonflux_nucleus_edff_nob(double& x, double& gamma, double& z, double& ra)
    __attribute__((weak));
}

// photon flux modelling switches (public module variables, default-kind LOGICALs)
extern "C" {
extern int __elasticphotonphotonflux_MOD_use_chargeformfactor4photon __attribute__((weak));
extern int __elasticphotonphotonflux_MOD_proton_ff_correction __attribute__((weak));
}

// two-photon fluxes (same names in all versions)
extern "C" {
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pp(double& x1, double& x2, int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pa_hardsphere(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_pa_woodssaxon(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_ab_hardsphere(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
extern double __elasticphotonphotonflux_MOD_photonphotonflux_ab_woodssaxon(double& x1,
                                                                           double& x2,
                                                                           int& force_p_no_hadronisation);
}

// gammaUPC version-dependent beams content -> resolve when loading
extern "C" {
extern gammaUPC::to_heavyion_pdg_t to_heavyion_pdg_ __attribute__((weak));
extern int __elasticphotonphotonflux_MOD_nucleara_beam1 __attribute__((weak));
extern int __elasticphotonphotonflux_MOD_nucleara_beam2 __attribute__((weak));
extern int __elasticphotonphotonflux_MOD_nuclearz_beam1 __attribute__((weak));
extern int __elasticphotonphotonflux_MOD_nuclearz_beam2 __attribute__((weak));
}

namespace gammaUPC {
  std::ostream& operator<<(std::ostream& os, const BeamsContent& content) {
    switch (content) {
      case BeamsContent::none:
      default:
        return os << "none (pp mode)";
      case BeamsContent::heavyion_cb:
        return os << "to_heavyion_pdg common block";
      case BeamsContent::module_var:
        return os << "ElasticPhotonPhotonFlux module variables";
    };
  }

  BeamsContent setBeamsContent(const std::array<int, 2>& z, const std::array<int, 2>& a) {
    if (&to_heavyion_pdg_) {  // (nb_proton, nb_neutron) = (1, 0) for protons
      for (size_t i = 0; i < 2; ++i)
        to_heavyion_pdg_.nb_proton[i] = z[i], to_heavyion_pdg_.nb_neutron[i] = a[i] - z[i];
      return BeamsContent::heavyion_cb;
    }
    if (&__elasticphotonphotonflux_MOD_nucleara_beam1 && &__elasticphotonphotonflux_MOD_nucleara_beam2 &&
        &__elasticphotonphotonflux_MOD_nuclearz_beam1 && &__elasticphotonphotonflux_MOD_nuclearz_beam2) {
      // convention in these versions: A = Z = 0 for a proton beam
      __elasticphotonphotonflux_MOD_nucleara_beam1 = a[0] > 1 ? a[0] : 0;
      __elasticphotonphotonflux_MOD_nuclearz_beam1 = a[0] > 1 ? z[0] : 0;
      __elasticphotonphotonflux_MOD_nucleara_beam2 = a[1] > 1 ? a[1] : 0;
      __elasticphotonphotonflux_MOD_nuclearz_beam2 = a[1] > 1 ? z[1] : 0;
      return BeamsContent::module_var;
    }
    if (a[0] > 1 || a[1] > 1)
      throw CG_FATAL("gammaUPC:setBeamsContent")
          << "Could not find any way to pass the ions content to this gammaUPC version (neither to_heavyion_pdg common "
             "block nor nuclearA_beamN/nuclearZ_beamN module variables).";
    return BeamsContent::none;
  }

  namespace {
    struct FormFactorSettings {
      bool set{false};  ///< was the modelling explicitly selected?
      FormFactor form_factor{FormFactor::EDFF};
      bool proton_ff_correction{false};
    } form_factor_settings;
  }  // namespace

  std::ostream& operator<<(std::ostream& os, const FormFactor& form_factor) {
    switch (form_factor) {
      case FormFactor::EDFF:
        return os << "EDFF";
      case FormFactor::ChFF:
        return os << "ChFF";
    }
    return os << "invalid{" << static_cast<int>(form_factor) << "}";
  }

  void setFormFactor(FormFactor form_factor, bool proton_ff_correction) {
    auto& settings = form_factor_settings;
    if (settings.set) {
      if (settings.form_factor != form_factor || settings.proton_ff_correction != proton_ff_correction)
        throw CG_FATAL("gammaUPC:setFormFactor")
            << "Photon flux modelling was already set to " << settings.form_factor
            << " (proton FF correction: " << settings.proton_ff_correction << "), and cannot be changed to "
            << form_factor << " (proton FF correction: " << proton_ff_correction
            << ") within the same process, as gammaUPC caches its fluxes grids.";
      return;
    }
    if (!&__elasticphotonphotonflux_MOD_use_chargeformfactor4photon) {
      if (form_factor != FormFactor::EDFF)
        throw CG_FATAL("gammaUPC:setFormFactor")
            << "This gammaUPC version does not allow to select the photon flux modelling (only EDFF available).";
    } else
      __elasticphotonphotonflux_MOD_use_chargeformfactor4photon = form_factor == FormFactor::ChFF ? 1 : 0;
    if (proton_ff_correction) {
      if (form_factor != FormFactor::ChFF)
        throw CG_FATAL("gammaUPC:setFormFactor") << "Proton form factor correction is only defined for ChFF.";
      if (!&__elasticphotonphotonflux_MOD_proton_ff_correction)
        throw CG_FATAL("gammaUPC:setFormFactor")
            << "This gammaUPC version does not allow to include the proton form factor correction.";
    }
    if (&__elasticphotonphotonflux_MOD_proton_ff_correction)
      __elasticphotonphotonflux_MOD_proton_ff_correction = proton_ff_correction ? 1 : 0;
    settings = FormFactorSettings{true, form_factor, proton_ff_correction};
    CG_DEBUG("gammaUPC:setFormFactor") << "Photon flux modelling set to " << form_factor
                                       << " (proton FF correction: " << proton_ff_correction << ").";
  }

  FormFactor formFactor() { return form_factor_settings.form_factor; }

  double photonFluxP(double x, double gamma) { return photonFluxP(x, gamma, form_factor_settings.form_factor); }

  double photonFluxP(double x, double gamma, FormFactor form_factor) {
    double flux{0.};
    switch (form_factor) {
      case FormFactor::EDFF: {
        static const auto flux_function = &__elasticphotonphotonflux_MOD_photonflux_proton_edff_nob
                                              ? &__elasticphotonphotonflux_MOD_photonflux_proton_edff_nob
                                              : &__elasticphotonphotonflux_MOD_photonflux_proton_nob;
        if (!flux_function)
          throw CG_FATAL("gammaUPC:photonFluxP") << "No single-photon EDFF proton flux found in this gammaUPC version.";
        flux = flux_function(x, gamma);
      } break;
      case FormFactor::ChFF: {
        if (!&__elasticphotonphotonflux_MOD_photonflux_proton_chff_nob)
          throw CG_FATAL("gammaUPC:photonFluxP") << "No single-photon ChFF proton flux found in this gammaUPC version.";
        int correction = form_factor_settings.proton_ff_correction ? 1 : 0;
        flux = __elasticphotonphotonflux_MOD_photonflux_proton_chff_nob(x, gamma, &correction);
      } break;
    }
    CG_DEBUG_LOOP("gammaUPC:photonFluxP") << "Single-photon flux from proton (" << form_factor << ", gamma=" << gamma
                                          << ", no b-effect): " << flux << ".";
    return flux;
  }

  double photonFluxA(double a, double gamma, double z, double ra) {
    if (form_factor_settings.form_factor != FormFactor::EDFF)
      throw CG_FATAL("gammaUPC:photonFluxA")
          << "Single-photon nucleus flux is only interfaced for the EDFF modelling (ChFF requires the Woods-Saxon "
             "parameters of the nucleus).";
    static const auto flux_function = &__elasticphotonphotonflux_MOD_photonflux_nucleus_edff_nob
                                          ? &__elasticphotonphotonflux_MOD_photonflux_nucleus_edff_nob
                                          : &__elasticphotonphotonflux_MOD_photonflux_nucleus_nob;
    if (!flux_function)
      throw CG_FATAL("gammaUPC:photonFluxA") << "No single-photon nucleus flux found in this gammaUPC version.";
    const auto flux = flux_function(a, gamma, z, ra);
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
            __elasticphotonphotonflux_MOD_photonphotonflux_pa_hardsphere(x1, x2, force_p_no_hadronisation);
        CG_DEBUG_LOOP("gammaUPC:twoPhotonFluxPA")
            << "Two-photon flux from proton-HI with hard sphere modelling for x1=" << x1 << ", x2=" << x2
            << " [non-had.: " << static_cast<bool>(force_p_no_hadronisation) << "]: " << flux << ".";
        return flux;
      }
      case HeavyIonMode::WoodsSaxon: {
        const auto flux =
            __elasticphotonphotonflux_MOD_photonphotonflux_pa_woodssaxon(x1, x2, force_p_no_hadronisation);
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
            __elasticphotonphotonflux_MOD_photonphotonflux_ab_hardsphere(x1, x2, force_p_no_hadronisation);
        CG_DEBUG_LOOP("gammaUPC:twoPhotonFluxAB")
            << "Two-photon flux from HI-HI with hard sphere modelling for x1=" << x1 << ", x2=" << x2
            << " [non-had.: " << static_cast<bool>(force_p_no_hadronisation) << "]: " << flux << ".";
        return flux;
      }
      case HeavyIonMode::WoodsSaxon: {
        const auto flux =
            __elasticphotonphotonflux_MOD_photonphotonflux_ab_woodssaxon(x1, x2, force_p_no_hadronisation);
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
