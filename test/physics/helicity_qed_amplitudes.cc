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

#include <array>
#include <cmath>
#include <vector>

#include "CepGen/Core/ParametersList.h"
#include "CepGen/HelicityAmplitudes/HelicityAmplitudes.h"
#include "CepGen/Utils/Test.h"

using namespace cepgen;

/// Check that the QED helicity amplitudes, once summed and averaged over helicities, reproduce the textbook
/// on-shell γγ → f f̄ and γγ → S⁺S⁻ squared matrix elements (as used in the on-shell methods of PPtoFF and PPtoSS)
int main() {
  const auto values = std::vector<std::array<double, 3> >{
      {100., 3., 0.3}, {100., 3., 0.9}, {100., 3., 2.1}, {1.e4, 0.1, 1.2}, {1.e4, 100., 0.7}};

  {  // fermions
    const auto fermion = amplitudes::HelicityAmplitudesFactory::get().build(ParametersList().setName("fermionQED"));
    CG_TEST(fermion->hasAAFF() && !fermion->hasAASS(), "fermion model capabilities");

    for (const auto& [s_hat, mass, theta] : values) {
      const auto m2 = mass * mass, beta = std::sqrt(1. - 4. * m2 / s_hat), e = 0.5 * std::sqrt(s_hat), p = e * beta,
                 cth = std::cos(theta);
      const auto t_hat = m2 - 2. * e * (e - p * cth), u_hat = m2 - 2. * e * (e + p * cth);
      const amplitudes::Kinematics kin(m2, s_hat, t_hat, u_hat);

      double sum_f{0.};
      for (const short l1 : {-1, 1})
        for (const short l2 : {-1, 1})
          for (const short s1 : {-1, 1})
            for (const short s2 : {-1, 1})
              sum_f += std::norm(fermion->amplitude(kin, {l1, l2, s1, s2}));
      const auto m4 = m2 * m2, m8 = m4 * m4;
      const auto out = 6. * m8 - 3. * m4 * t_hat * t_hat - 14. * m4 * t_hat * u_hat - 3. * m4 * u_hat * u_hat +
                       m2 * t_hat * t_hat * t_hat + 7. * m2 * t_hat * t_hat * u_hat + 7. * m2 * t_hat * u_hat * u_hat +
                       m2 * u_hat * u_hat * u_hat - t_hat * t_hat * t_hat * u_hat - t_hat * u_hat * u_hat * u_hat;
      const auto ref_f = -2. * out / std::pow((m2 - t_hat) * (m2 - u_hat), 2);  // = <|M|^2> for e = Q = 1
      CG_TEST_EQUIV(0.25 * sum_f / ref_f, 1., "gamma-gamma -> f-fbar helicity sum");
    }

    {  // scalars
      const auto scalar = amplitudes::HelicityAmplitudesFactory::get().build(ParametersList().setName("scalarQED"));
      CG_TEST(scalar->hasAASS() && !scalar->hasAAFF(), "scalar model capabilities");

      for (const auto& [s_hat, mass, theta] : values) {
        const auto m2 = mass * mass, beta = std::sqrt(1. - 4. * m2 / s_hat), e = 0.5 * std::sqrt(s_hat), p = e * beta,
                   cth = std::cos(theta);
        const auto t_hat = m2 - 2. * e * (e - p * cth), u_hat = m2 - 2. * e * (e + p * cth);
        const amplitudes::Kinematics kin(m2, s_hat, t_hat, u_hat);
        double sum_s{0.};
        for (const short l1 : {-1, 1})
          for (const short l2 : {-1, 1})
            sum_s += std::norm(scalar->amplitude(kin, {l1, l2, 0, 0}));
        const auto d = 4. * (t_hat - m2) * (u_hat - m2) / (s_hat * s_hat), omb2 = 4. * m2 / s_hat, b2s2 = d - omb2;
        const auto ref_s = 2. * (omb2 * omb2 + b2s2 * b2s2) / (d * d);
        CG_TEST_EQUIV(0.25 * sum_s / ref_s, 1., "gamma-gamma -> s-sbar helicity sum");
      }
    }
  }
  CG_TEST_SUMMARY;
}
