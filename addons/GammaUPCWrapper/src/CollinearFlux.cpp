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
#include <CepGen/Physics/Beam.h>
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
          q2_range_(steer<Limits>("q2range")),
          inv_log_q2_range_(1. / std::log(q2_range_.max() / q2_range_.min())) {
      if (q2_range_.min() <= 0. || q2_range_.max() <= q2_range_.min())
        throw CG_FATAL("gammaUPC:CollinearFlux") << "Invalid virtuality range: " << q2_range_ << ".";
    }

    static ParametersDescription description() {
      auto desc = cepgen::CollinearFlux::description();
      desc.setDescription("gammaUPC single-beam photon flux");
      desc.add("q2range", Limits{1.e-10, 5.})
          .setDescription("parton virtuality range sampled by the collinear parton sampler (GeV^2)");
      return desc;
    }

    void initialise(const Beam& beam) override {  // beam particle content and Lorentz factor from beam kinematics
      const auto& momentum = beam.momentum();
      if (mass2_ = momentum.mass2(); mass2_ <= 0. || (gamma_ = momentum.energy() / momentum.mass()) <= 1.)
        throw CG_FATAL("gammaUPC:CollinearFlux") << "Invalid beam kinematics: " << momentum << ".";
      if (heavy_ion_ = HeavyIon::isHI(beam.integerPdgId()); heavy_ion_) {
        const auto hi = HeavyIon::fromPdgId(beam.integerPdgId());
        z_ = static_cast<double>(hi.Z);
        a_ = static_cast<double>(hi.A);
        radius_fm_ = hi.radius() * 1.e15;  // m -> fm
      }
      CG_DEBUG("gammaUPC:CollinearFlux") << "Initialised for beam " << beam.integerPdgId() << " with gamma=" << gamma_
                                         << (heavy_ion_ ? " (heavy ion)" : "") << ".";
    }

    bool fragmenting() const override { return false; }
    spdgid_t partonPdgId() const override { return PDG::photon; }
    double mass2() const override { return mass2_; }

    double fluxQ2(double x, double q2) const override {
      if (gamma_ <= 0.)
        throw CG_FATAL("gammaUPC:CollinearFlux") << "Flux was not initialised with a beam.";
      if (!x_range_.contains(x, true) || !q2_range_.contains(q2))
        return 0.;
      if (heavy_ion_)
        return a_ * x * photonFluxA(x, gamma_, z_, radius_fm_) * inv_log_q2_range_;
      return x * photonFluxP(x, gamma_) * inv_log_q2_range_;
    }

  private:
    const Limits q2_range_;
    const double inv_log_q2_range_;

    bool heavy_ion_{false};
    double gamma_{0.}, z_{1.}, a_{1.}, radius_fm_{0.}, mass2_{0.};
  };
}  // namespace gammaUPC
using gammaUPCCollinearFlux = gammaUPC::CollinearFlux;
REGISTER_COLLINEAR_FLUX("gammaUPC", gammaUPCCollinearFlux);
