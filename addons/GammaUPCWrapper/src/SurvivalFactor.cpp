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
          heavy_ion_mode_(steerAs<int, HeavyIonMode>("heavyIonMode")) {}

    static ParametersDescription description() {
      auto desc = cepgen::SurvivalFactor::description();
      desc.setDescription("gammaUPC two-photon fluxes non-overlap correction");
      desc.add("forcePNonHadronisation", true)
          .setDescription("disable the hadronisation of beam remnants after parton emission?");
      desc.addAs<int, HeavyIonMode>("heavyIonMode", HeavyIonMode::HardSphere)
          .setDescription("type of modelling for the heavy ion parton emission")
          .allow(static_cast<int>(HeavyIonMode::HardSphere), "hard sphere model")
          .allow(static_cast<int>(HeavyIonMode::WoodsSaxon), "Woods-Saxon model");
      return desc;
    }

    double operator()(const Event& event) const override {
      if (!initialised_)
        initialise(event);
      // parton momentum fractions, as computed by the factorised process from the central system
      double x1 = 0., x2 = 0.;
      for (const auto& part : event(Particle::Role::CentralSystem)) {
        const auto& mom = part.momentum();
        x1 += std::fabs(mom.energy() + mom.pz());
        x2 += std::fabs(mom.energy() - mom.pz());
      }
      x1 /= sqrt_s_, x2 /= sqrt_s_;
      if (x1 <= 0. || x1 >= 1. || x2 <= 0. || x2 >= 1.)
        return 0.;
      if (const auto single_fluxes = singleFlux(0, x1) * singleFlux(1, x2); single_fluxes > 0.)
        return twoPhotonFlux(x1, x2) / single_fluxes;
      return 0.;
    }

  private:
    enum struct BeamsMode { pp, pB, Ap, AB };

    void initialise(const Event& event) const {
      const std::array<const Particle*, 2> beams{&event.oneWithRole(Particle::Role::IncomingBeam1),
                                                 &event.oneWithRole(Particle::Role::IncomingBeam2)};
      std::array<int, 2> z{1, 1}, a{1, 1};
      for (size_t i = 0; i < 2; ++i) {
        const auto& mom = beams.at(i)->momentum();
        gamma_.at(i) = mom.energy() / mom.mass();
        if (heavy_ion_.at(i) = HeavyIon::isHI(beams.at(i)->integerPdgId()); heavy_ion_.at(i)) {
          const auto hi = HeavyIon::fromPdgId(beams.at(i)->integerPdgId());
          z.at(i) = static_cast<int>(hi.Z), a.at(i) = static_cast<int>(hi.A);
          radius_fm_.at(i) = hi.radius() * 1.e15;  // m -> fm
        }
      }
      if (heavy_ion_.at(0))
        beams_mode_ = heavy_ion_.at(1) ? BeamsMode::AB : BeamsMode::Ap;
      else
        beams_mode_ = heavy_ion_.at(1) ? BeamsMode::pB : BeamsMode::pp;
      z_ = {static_cast<double>(z.at(0)), static_cast<double>(z.at(1))};
      if (beams_mode_ == BeamsMode::Ap) {  // gammaUPC expects proton as first beam for asymmetric collisions
        std::swap(z[0], z[1]);
        std::swap(a[0], a[1]);
      }
      const auto beams_content = setBeamsContent(z, a);
      to_collider_.ebeam[0] = beams[0]->momentum().energy();
      to_collider_.ebeam[1] = beams[1]->momentum().energy();
      if (beams_mode_ == BeamsMode::Ap)
        std::swap(to_collider_.ebeam[0], to_collider_.ebeam[1]);
      sqrt_s_ = event.cmEnergy();
      initialised_ = true;
      CG_DEBUG("gammaUPC:SurvivalFactor") << "Initialised with beams content passed through: " << beams_content
                                          << ", Lorentz factors: " << gamma_[0] << "/" << gamma_[1] << ".";
    }

    double singleFlux(size_t i, double x) const {
      if (heavy_ion_.at(i))
        return photonFluxA(x, gamma_.at(i), z_.at(i), radius_fm_.at(i));
      return photonFluxP(x, gamma_.at(i));
    }

    double twoPhotonFlux(double x1, double x2) const {
      switch (beams_mode_) {
        case BeamsMode::pp:
          return twoPhotonFluxPP(x1, x2, force_p_non_hadronisation_);
        case BeamsMode::pB:
          return twoPhotonFluxPA(x1, x2, force_p_non_hadronisation_, heavy_ion_mode_);
        case BeamsMode::Ap:
          return twoPhotonFluxPA(x2, x1, force_p_non_hadronisation_, heavy_ion_mode_);
        case BeamsMode::AB:
          return twoPhotonFluxAB(x1, x2, force_p_non_hadronisation_, heavy_ion_mode_);
      }
      throw CG_FATAL("gammaUPC:SurvivalFactor") << "Unknown beams mode: " << static_cast<int>(beams_mode_) << ".";
    }

    const bool force_p_non_hadronisation_;
    const HeavyIonMode heavy_ion_mode_;

    // beams information, retrieved from the first event
    mutable bool initialised_{false};
    mutable BeamsMode beams_mode_{BeamsMode::pp};
    mutable std::array<bool, 2> heavy_ion_{false, false};
    mutable std::array<double, 2> gamma_{0., 0.}, z_{1., 1.}, radius_fm_{0., 0.};
    mutable double sqrt_s_{0.};
  };
}  // namespace gammaUPC
using gammaUPCSurvivalFactor = gammaUPC::SurvivalFactor;
REGISTER_SURVIVAL_FACTOR("gammaUPC", gammaUPCSurvivalFactor);
