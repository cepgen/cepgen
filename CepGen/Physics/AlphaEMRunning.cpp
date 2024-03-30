/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022  Laurent Forthomme
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

#include <gsl/gsl_sf_zeta.h>

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/Message.h"

namespace cepgen {
  /// Running electromagnetic alpha calculator
  class AlphaEMRunning final : public Coupling {
  public:
    explicit AlphaEMRunning(const ParametersList& params) : Coupling(params), c_(steer<double>("c")) {}

    static ParametersDescription description() {
      auto desc = Coupling::description();
      desc.setDescription("Running alpha(EM) evolution algorithm");
      desc.add<double>("c", 1.).setDescription(
          "running parameter (0 is constant alphaQED, 1 is QED evolution, best L3 fit value is 1.05 +- 0.07 +- 0.14)");
      return desc;
    }

    double operator()(double q) const override {
      return constants::ALPHA_EM / (1. - c_ * (deltaAlphaL(q) + deltaAlphaH(q)));
    }

    double deltaAlphaL(double q) const {
      // lepton contribution to alphaEM (photon propagator), as parameterised by
      //   Steinhauser, Phys.Lett.B 429 (1998) 158-161
      //   https://doi.org/10.1016/S0370-2693(98)00503-6

      // only retrieve once lepton masses
      static const auto m2_el = std::pow(PDG::get().mass(11), 2), m2_mu = std::pow(PDG::get().mass(13), 2),
                        m2_tau = std::pow(PDG::get().mass(15), 2);
      // only compute once Riemann zeta function integer values
      static const auto zeta2 = gsl_sf_zeta_int(2), zeta3 = gsl_sf_zeta_int(3), zeta5 = gsl_sf_zeta_int(5);

      const auto q2 = q * q;

      //----- definition of all polarisation functions

      auto log_qm = [](double q2, double ml2) { return log(q2 / ml2); };
      // one-loop corrections
      auto pi0 = [&log_qm](double q2, double ml2) { return 20. / 9 - 4. / 3 * log_qm(q2, ml2) + 8. * ml2 / q2; };
      // two-loop corrections
      auto pi1 = [&log_qm](double q2, double ml2) {
        const auto lqm = log_qm(q2, ml2);
        return 5. / 6 - 4. * zeta3 - lqm - 12. * lqm * ml2 / q2;
      };
      // three-loop corrections
      auto pi2A = [&log_qm](double q2, double ml2) {  // quenched
        return -121. / 48 + (-5. + 8. * log(2)) * zeta2 - 99. / 16 * zeta3 + 10. * zeta5 + 0.125 * log_qm(q2, ml2);
      };
      auto pi2l = [&log_qm](double q2, double ml12, double ml22) {
        const auto lqm1 = log_qm(q2, ml12), lqm2 = log_qm(q2, ml22);
        return -116. / 27 + 4. / 3 * zeta2 + 38. / 9 * zeta3 + 14. / 9 * lqm1 + (5. / 18 - 4. / 3 * zeta3) * lqm2 +
               1. / 6 * lqm1 * lqm1 - lqm1 * lqm2 / 3.;
      };
      auto pi2F = [&log_qm](double q2, double ml2) {
        const auto lqm = log_qm(q2, ml2);
        return -307. / 216 - 8. / 3 * zeta2 + 545. / 144 * zeta3 + (11. / 6 - 4. / 3 * zeta3) * lqm - lqm * lqm / 6.;
      };
      auto pi2h = [&log_qm](double q2, double ml2) {
        const auto lqm = log_qm(q2, ml2);
        return -37. / 6 + 38. * zeta3 / 9. + (11. / 6 - 4. * zeta3 / 3.) * lqm - lqm * lqm / 6;
      };

      //----- compute lepton alphaQED contributions for all orders

      const auto alpha_ov_pi = constants::ALPHA_EM * M_1_PI;
      const auto order0 = -0.25 * alpha_ov_pi * (pi0(q2, m2_el) + pi0(q2, m2_mu) + pi0(q2, m2_tau));
      const auto order1 = -0.25 * alpha_ov_pi * alpha_ov_pi * (pi1(q2, m2_el) + pi1(q2, m2_mu) + pi1(q2, m2_tau));
      const auto order2_quenched = pi2A(q2, m2_el) + pi2A(q2, m2_mu) + pi2A(q2, m2_tau);
      const auto order2_l = pi2l(q2, m2_mu, m2_el) + pi2l(q2, m2_tau, m2_mu) + pi2l(q2, m2_tau, m2_el);
      const auto order2_F = pi2F(q2, m2_el) + pi2F(q2, m2_mu) + pi2F(q2, m2_tau);
      const auto order2_h = pi2h(q2, m2_el) + pi2h(q2, m2_mu) + pi2h(q2, m2_tau);
      const auto order2 =
          -0.25 * alpha_ov_pi * alpha_ov_pi * alpha_ov_pi * (order2_quenched + order2_l + order2_F + order2_h);

      return order0 + order1 + order2;
    }

    double deltaAlphaH(double q) const {
      // hadronic contribution to alphaEM, as parameterised by
      //   Burkhardt and Pietrzyk, Phys.Lett.B 513 (2001) 46-52
      //   https://doi.org/10.1016/S0370-2693(01)00393-8
      if (q < 0.)
        return 0.;
      auto param = [](double q2, double a, double b, double c) { return a + b * log1p(c * q2); };
      if (q <= 0.7)
        return param(q * q, 0., 0.0023092, 3.9925370);
      if (q <= 2.)
        return param(q * q, 0., 0.0022333, 4.2191779);
      if (q <= 4.)
        return param(q * q, 0., 0.0024402, 3.2496684);
      if (q <= 10.)
        return param(q * q, 0., 0.0027340, 2.0995092);
      static const double mZ = 91.1876;
      if (q <= mZ)
        return param(q * q, 0.0010485, 0.0029431, 1.);
      if (q <= 1.e4)
        return param(q * q, 0.0012234, 0.0029237, 1.);
      if (q <= 1.e5)
        return param(q * q, 0.0016894, 0.0028984, 1.);
      CG_WARNING("AlphaEMRunning:deltaAlpha") << "Q exceeds the validity range of Burkhardt et al. parameterisation.";
      return deltaAlphaH(1.e5);
    }

  private:
    const double c_;
  };
}  // namespace cepgen
REGISTER_ALPHAEM_MODULE("running", AlphaEMRunning);
