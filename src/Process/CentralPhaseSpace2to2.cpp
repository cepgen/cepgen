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

#include <cmath>

#include "CepGen/Modules/CentralPhaseSpaceFactory.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Process/CentralPhaseSpace.h"

using namespace cepgen;

/// Two-parton to two-particle central system phase space mapper
/// \note Variables: both central particles' rapidities, and their transverse momentum difference
///  (norm and azimuthal angle)
class CentralPhaseSpace2to2 final : public CentralPhaseSpace {
public:
  explicit CentralPhaseSpace2to2(const ParametersList& params) : CentralPhaseSpace(params) {}

  static ParametersDescription description();

  size_t multiplicity() const override { return 2; }
  std::vector<IntegrationVariable> integrationVariables() override;
  double generate(const Momentum& qt_sum, std::vector<Momentum>& central_momenta) const override;

private:
  // factor 1/4 from jacobian of transformations
  static constexpr double prefactor_ = M_1_PI * M_1_PI / 16. / 4.;
  // mapped variables
  double m_y_c1_{0.};         ///< Rapidity of the first central particle
  double m_y_c2_{0.};         ///< Rapidity of the second central particle
  double m_pt_diff_{0.};      ///< Transverse momentum difference for the two central particle
  double m_phi_pt_diff_{0.};  ///< Azimuthal angle difference for the two central particles
};

ParametersDescription CentralPhaseSpace2to2::description() {
  auto desc = CentralPhaseSpace::description();
  desc.setDescription("2-to-2 central system phase space mapper");
  return desc;
}

std::vector<IntegrationVariable> CentralPhaseSpace2to2::integrationVariables() {
  const auto& kin_cuts = kinematics().cuts().central;
  const auto lim_rap = kin_cuts.rapidity_single.truncate(Limits{-6., 6.});
  return {{m_y_c1_, proc::Process::Mapping::linear, lim_rap, "y1", "First outgoing particle rapidity"},
          {m_y_c2_, proc::Process::Mapping::linear, lim_rap, "y2", "Second outgoing particle rapidity"},
          {m_pt_diff_,
           proc::Process::Mapping::linear,
           kin_cuts.pt_diff.truncate(Limits{0., 500.}),
           "pt_diff",
           "Final state particles transverse momentum difference"},
          {m_phi_pt_diff_,
           proc::Process::Mapping::linear,
           kin_cuts.phi_diff.truncate(Limits{0., 2. * M_PI}),
           "phi_pt_diff",
           "Final state particles azimuthal angle difference"}};
}

double CentralPhaseSpace2to2::generate(const Momentum& qt_sum, std::vector<Momentum>& central_momenta) const {
  // central system cuts are read at each point, as processes may still modify them after initialisation
  const auto& cuts = kinematics().cuts().central;
  if (!cuts.rapidity_diff.contains(std::fabs(m_y_c1_ - m_y_c2_)))  // rapidity distance
    return 0.;
  const auto pt_diff = Momentum::fromPtEtaPhiE(m_pt_diff_, 0., m_phi_pt_diff_);
  const auto pt_c1 = 0.5 * (qt_sum + pt_diff), pt_c2 = 0.5 * (qt_sum - pt_diff);
  const auto p1t = pt_c1.pt(), p2t = pt_c2.pt();
  // apply user cuts on central system
  if (!cuts.pt_single.contains(p1t) || !particleCuts(0).pt_single.contains(p1t))
    return 0.;
  if (!cuts.pt_single.contains(p2t) || !particleCuts(1).pt_single.contains(p2t))
    return 0.;
  if (!cuts.pt_diff.contains(std::fabs(p1t - p2t)))  // transverse momentum difference
    return 0.;
  const auto phi1 = pt_c1.phi(), phi2 = pt_c2.phi();
  if (!cuts.phi_single.contains(phi1) || !particleCuts(0).phi_single.contains(phi1))
    return 0.;
  if (!cuts.phi_single.contains(phi2) || !particleCuts(1).phi_single.contains(phi2))
    return 0.;
  // four-momenta of the outgoing central particles
  central_momenta.resize(2);
  central_momenta[0] = Momentum::fromPtYPhiM(p1t, m_y_c1_, phi1, mass(0));
  central_momenta[1] = Momentum::fromPtYPhiM(p2t, m_y_c2_, phi2, mass(1));
  return prefactor_ * m_pt_diff_;
}
REGISTER_CENTRAL_PHASE_SPACE("2to2", CentralPhaseSpace2to2);
