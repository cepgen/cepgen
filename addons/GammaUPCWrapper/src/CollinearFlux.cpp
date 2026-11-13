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
#include <CepGen/Modules/PartonFluxFactory.h>
#include <CepGen/PartonFluxes/CollinearFlux.h>
#include <CepGen/Physics/HeavyIon.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/Limits.h>

#include <cmath>

#include "CepGenGammaUPC/GammaUPCInterface.h"

using namespace cepgen;

namespace gammaUPC {
  /// Single-beam gammaUPC photon flux (EDFF, no impact parameter dependence)
  class CollinearFlux final : public cepgen::CollinearFlux {
  public:
    explicit CollinearFlux(const ParametersList& params)
        : cepgen::CollinearFlux(params),
          beam_pdg_id_(steer<pdgid_t>("pdgId")),
          heavy_ion_(HeavyIon::isHI(beam_pdg_id_)),
          gamma_(steer<double>("gamma")),
          q2_range_(steer<Limits>("q2range")),
          inv_log_q2_range_(1. / std::log(q2_range_.max() / q2_range_.min())) {
      if (gamma_ <= 1.)
        throw CG_FATAL("gammaUPC:CollinearFlux") << "Invalid beam Lorentz factor: " << gamma_ << ".";
      if (q2_range_.min() <= 0. || q2_range_.max() <= q2_range_.min())
        throw CG_FATAL("gammaUPC:CollinearFlux") << "Invalid virtuality range: " << q2_range_ << ".";
      if (heavy_ion_) {
        const auto hi = HeavyIon::fromPdgId(beam_pdg_id_);
        z_ = static_cast<double>(hi.Z);
        a_ = static_cast<double>(hi.A);
        radius_fm_ = hi.radius() * 1.e15;  // m -> fm
        mass2_ = std::pow(hi.mass(), 2);
      } else
        mass2_ = std::pow(PDG::get().mass(beam_pdg_id_), 2);
    }

    static ParametersDescription description() {
      auto desc = cepgen::CollinearFlux::description();
      desc.setDescription("gammaUPC single-beam photon flux");
      desc.addAs<int, pdgid_t>("pdgId", PDG::proton).setDescription("beam particle (proton or heavy ion) PDG id");
      desc.add("gamma", 0.).setDescription("beam Lorentz factor");
      desc.add("q2range", Limits{1.e-10, 5.})
          .setDescription("parton virtuality range sampled by the collinear parton sampler (GeV^2)");
      return desc;
    }

    bool fragmenting() const override { return false; }
    spdgid_t partonPdgId() const override { return PDG::photon; }
    double mass2() const override { return mass2_; }

    double fluxQ2(double x, double q2) const override {
      if (!x_range_.contains(x, true) || !q2_range_.contains(q2))
        return 0.;
      if (heavy_ion_)
        return a_ * x * photonFluxA(x, gamma_, z_, radius_fm_) * inv_log_q2_range_;
      return x * photonFluxP(x, gamma_) * inv_log_q2_range_;
    }

  private:
    const pdgid_t beam_pdg_id_;
    const bool heavy_ion_;
    const double gamma_;
    const Limits q2_range_;
    const double inv_log_q2_range_;

    double z_{1.}, a_{1.}, radius_fm_{0.}, mass2_{0.};
  };
}  // namespace gammaUPC
using gammaUPCCollinearFlux = gammaUPC::CollinearFlux;
REGISTER_COLLINEAR_FLUX("gammaUPC", gammaUPCCollinearFlux);
