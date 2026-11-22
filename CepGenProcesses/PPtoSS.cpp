/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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
#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/FactorisedProcess.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/RandomGenerator.h"

using namespace cepgen;
using namespace std::string_literals;

/// Compute the matrix element for a CE \f$\gamma\gamma\rightarrow S^+S^-\f$ process, with \f$S\f$ a charged scalar
class PPtoSS final : public proc::FactorisedProcess {
public:
  explicit PPtoSS(const ParametersList& params)
      : proc::FactorisedProcess(params, make_pdgids_pair(params.get<ParticleProperties>("pair").pdgid)) {}

  proc::ProcessPtr clone() const override { return std::make_unique<PPtoSS>(parameters()); }

  static ParametersDescription description() {
    auto desc = FactorisedProcess::description();
    desc.setDescription("γγ → S⁺S⁻ (charged scalars pair, e.g. sleptons)");
    desc.addAs<int, pdgid_t>("pair", PDG::piPlus)
        .setDescription("type of central charged scalar emitted (must be defined in the PDG library)");
    return desc;
  }

private:
  static spdgids_t make_pdgids_pair(pdgid_t pair) {
    return spdgids_t{static_cast<spdgid_t>(pair), -static_cast<spdgid_t>(pair)};
  }

  void prepareFactorisedPhaseSpace() override {
    const auto& cs_prop = PDG::get()(central().at(0));
    if (cs_prop.fermion)
      throw CG_FATAL("PPtoSS:prepare") << "Invalid central particle selected (fermion; use 'pptoff' instead): "
                                       << cs_prop << ".";
    if (cs_prop.integerCharge() == 0)
      throw CG_FATAL("PPtoSS:prepare") << "The central scalar particle must be charged: " << cs_prop << ".";
    if (cs_prop.mass <= 0.)
      throw CG_FATAL("PPtoSS:prepare") << "The central scalar particle must be massive: " << cs_prop << ".";
    for (const auto role : {Particle::Role::Parton1, Particle::Role::Parton2})
      if (const auto parton_id = event().oneWithRole(role).pdgId(); parton_id != PDG::photon)
        throw CG_FATAL("PPtoSS:prepare") << "Unsupported parton id: '" << parton_id << "'. Only photons are allowed.";

    ms2_ = cs_prop.mass * cs_prop.mass;
    qs2_ = static_cast<double>(cs_prop.integerCharge() * cs_prop.integerCharge()) / 9.;
    // sum over the final state colours (scalar leptons-like particles have no colour)
    colour_factor_ = cs_prop.colours > 0 ? cs_prop.colours : 1.;

    CG_DEBUG("PPtoSS:prepare") << "Incoming beams: mA = " << mA() << " GeV/mB = " << mB() << " GeV.\n\t"
                               << "Produced particles: " << central() << ".";

    // constrain central particles cuts
    if (!kinematics().cuts().central.pt_diff.valid())
      kinematics().cuts().central.pt_diff = {0., 50.};
  }

  double computeFactorisedMatrixElement() override {
    const double s_hat = shat(), t_hat = that(), u_hat = uhat();            // buffer Mandelstam variables
    const auto d = 4. * (t_hat - ms2_) * (u_hat - ms2_) / (s_hat * s_hat);  // = 1-beta^2 cos^2(theta)
    if (!utils::positive(d))
      return 0.;
    const auto q = std::sqrt(std::fabs(t_hat));  // t-hat negative in region of interest
    const auto prefactor = colour_factor_ * kFourPi * kFourPi * qs2_ * qs2_ * alphaEM(q) * alphaEM(q);
    if (!utils::positive(prefactor))
      return 0.;

    const auto one_minus_beta2 = 4. * ms2_ / s_hat, beta2_sin2 = d - one_minus_beta2;
    return 2. * prefactor * (one_minus_beta2 * one_minus_beta2 + beta2_sin2 * beta2_sin2) / (d * d) / (s_hat * s_hat);
  }

  static constexpr double kFourPi = 4. * M_PI;
  double ms2_{0.}, qs2_{0.}, colour_factor_{1.};
};
REGISTER_PROCESS("pptoss", PPtoSS);
