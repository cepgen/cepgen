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
#include <CepGen/Event/Event.h>
#include <CepGen/Modules/SurvivalFactorFactory.h>
#include <CepGen/Physics/HeavyIon.h>
#include <CepGen/SurvivalFactors/SurvivalFactor.h>

#include <cmath>

#include "CepGenGammaUPC/GammaUPCInterface.h"

using namespace cepgen;

namespace gammaUPC {
  /// Beam-beam non-overlap correction from the gammaUPC two-photon fluxes
  class SurvivalFactor final : public cepgen::SurvivalFactor {
  public:
    explicit SurvivalFactor(const ParametersList& params)
        : cepgen::SurvivalFactor(params),
          force_p_non_hadronisation_(steer<bool>("forcePNonHadronisation")),
          heavy_ion_mode_(steerAs<int, HeavyIonMode>("heavyIonMode")),
          form_factor_(steerAs<int, FormFactor>("formFactor")),
          proton_ff_correction_(steer<bool>("protonFFCorrection")) {}

    static ParametersDescription description() {
      auto desc = cepgen::SurvivalFactor::description();
      desc.setDescription("gammaUPC two-photon fluxes non-overlap correction");
      desc.add("forcePNonHadronisation", false)
          .setDescription(
              "force the no-hadronic-interaction probability to unity (i.e. switch off the survival factor)?");
      desc.addAs<int, HeavyIonMode>("heavyIonMode", HeavyIonMode::HardSphere)
          .setDescription("type of modelling for the heavy ion parton emission")
          .allow(static_cast<int>(HeavyIonMode::HardSphere), "hard sphere model")
          .allow(static_cast<int>(HeavyIonMode::WoodsSaxon), "Woods-Saxon model");
      desc.addAs<int, FormFactor>("formFactor", FormFactor::EDFF)
          .setDescription("photon flux modelling (process-wide gammaUPC setting)")
          .allow(static_cast<int>(FormFactor::EDFF), "electric dipole form factor")
          .allow(static_cast<int>(FormFactor::ChFF), "charge form factor (Woods-Saxon heavy ion modelling only)");
      desc.add("protonFFCorrection", false)
          .setDescription("include the proton form factor correction beyond the dipole (ChFF only)?");
      return desc;
    }

    double operator()(const Event& event) const override {
      if (!initialised_)
        initialise(event);
      double q_plus = 0., q_minus = 0.;
      for (const auto& part : event(Particle::Role::CentralSystem)) {
        const auto& mom = part.momentum();
        q_plus += std::fabs(mom.energy() + mom.pz());
        q_minus += std::fabs(mom.energy() - mom.pz());
      }
      const auto x1 = mass_numbers_.at(0) * q_plus / beams_light_cone_.at(0),
                 x2 = mass_numbers_.at(1) * q_minus / beams_light_cone_.at(1);
      if (x1 <= 0. || x1 >= 1. || x2 <= 0. || x2 >= 1.)
        return 0.;
      if (const auto no_overlap_flux = twoPhotonFlux(x1, x2, true); no_overlap_flux > 0.)
        return twoPhotonFlux(x1, x2, force_p_non_hadronisation_) / no_overlap_flux;
      return 0.;
    }

  private:
    enum struct BeamsMode { pp, pB, Ap, AB };

    void initialise(const Event& event) const {
      const std::array<const Particle*, 2> beams{&event.oneWithRole(Particle::Role::IncomingBeam1),
                                                 &event.oneWithRole(Particle::Role::IncomingBeam2)};
      std::array<int, 2> z{1, 1}, a{1, 1};
      std::array<bool, 2> heavy_ion{false, false};
      for (size_t i = 0; i < 2; ++i) {
        const auto& mom = beams.at(i)->momentum();
        beams_light_cone_.at(i) = mom.energy() + (i == 0 ? +1. : -1.) * mom.pz();
        if (heavy_ion.at(i) = HeavyIon::isHI(beams.at(i)->integerPdgId()); heavy_ion.at(i)) {
          const auto hi = HeavyIon::fromPdgId(beams.at(i)->integerPdgId());
          z.at(i) = static_cast<int>(hi.Z), a.at(i) = static_cast<int>(hi.A);
        }
        mass_numbers_.at(i) = a.at(i);
      }
      if (heavy_ion.at(0))
        beams_mode_ = heavy_ion.at(1) ? BeamsMode::AB : BeamsMode::Ap;
      else
        beams_mode_ = heavy_ion.at(1) ? BeamsMode::pB : BeamsMode::pp;
      if (form_factor_ == FormFactor::ChFF && beams_mode_ != BeamsMode::pp &&
          heavy_ion_mode_ != HeavyIonMode::WoodsSaxon)
        throw CG_FATAL("gammaUPC:SurvivalFactor")
            << "ChFF photon fluxes are only available with the Woods-Saxon heavy ion modelling in gammaUPC.";
      setFormFactor(form_factor_, proton_ff_correction_);  // before any flux evaluation (cached grids)
      if (beams_mode_ == BeamsMode::Ap) {  // gammaUPC expects proton as first beam for asymmetric collisions
        std::swap(z[0], z[1]);
        std::swap(a[0], a[1]);
      }
      const auto beams_content = setBeamsContent(z, a);
      to_collider_.ebeam[0] = beams[0]->momentum().energy();
      to_collider_.ebeam[1] = beams[1]->momentum().energy();
      if (beams_mode_ == BeamsMode::Ap)
        std::swap(to_collider_.ebeam[0], to_collider_.ebeam[1]);
      initialised_ = true;
      CG_DEBUG("gammaUPC:SurvivalFactor")
          << "Initialised with beams content passed through: " << beams_content
          << ", beams energies: " << to_collider_.ebeam[0] << "/" << to_collider_.ebeam[1] << " GeV.";
    }

    double twoPhotonFlux(double x1, double x2, bool no_hadronic_interaction) const {
      const int no_had = no_hadronic_interaction ? 1 : 0;
      switch (beams_mode_) {
        case BeamsMode::pp:
          return twoPhotonFluxPP(x1, x2, no_had);
        case BeamsMode::pB:
          return twoPhotonFluxPA(x1, x2, no_had, heavy_ion_mode_);
        case BeamsMode::Ap:
          return twoPhotonFluxPA(x2, x1, no_had, heavy_ion_mode_);
        case BeamsMode::AB:
          return twoPhotonFluxAB(x1, x2, no_had, heavy_ion_mode_);
      }
      throw CG_FATAL("gammaUPC:SurvivalFactor") << "Unknown beams mode: " << static_cast<int>(beams_mode_) << ".";
    }

    const bool force_p_non_hadronisation_;
    const HeavyIonMode heavy_ion_mode_;
    const FormFactor form_factor_;
    const bool proton_ff_correction_;

    // beams information, retrieved from the first event
    mutable bool initialised_{false};
    mutable BeamsMode beams_mode_{BeamsMode::pp};
    mutable std::array<double, 2> beams_light_cone_{0., 0.};  ///< beams p^+ (positive-z) and p^- (negative-z)
    mutable std::array<double, 2> mass_numbers_{1., 1.};      ///< beams mass numbers (1 for protons)
  };
}  // namespace gammaUPC
using gammaUPCSurvivalFactor = gammaUPC::SurvivalFactor;
REGISTER_SURVIVAL_FACTOR("gammaUPC", gammaUPCSurvivalFactor);
