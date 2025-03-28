/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2021-2025  Laurent Forthomme
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

#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/Coupling.h"

using namespace cepgen;

/// Electromagnetic alpha running calculator
/// \note Shamelessly stolen from JETSET/PYTHIA
class AlphaEMBurkhardt final : public Coupling {
public:
  explicit AlphaEMBurkhardt(const ParametersList& params) : Coupling(params), q2min_(steer<double>("q2min")) {}

  static ParametersDescription description() {
    auto desc = Coupling::description();
    desc.setDescription("Burkhardt et al. alpha(EM) evolution algorithm");
    desc.add("q2min", 2.e-6).setDescription("Minimum Q^2 to start alpha(EM) evolution");
    return desc;
  }

  double operator()(double q) const override {
    const double q2 = q * q;
    if (q2 < q2min_)
      return constants::ALPHA_EM;
    const double log_q2 = log(q2), log_1_pl_q2 = log1p(q2);
    // Calculate real part of photon vacuum polarization.
    // - for leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
    // - for hadrons use parametrization of H. Burkhardt et al.
    // See R. Kleiss et al., CERN 89-08, vol. 3, pp. 129-131.
    double rpigg;
    if (q2 < 9.e-2)
      rpigg = AEM_3PI * (13.4916 + log_q2) + 0.00835 * log_1_pl_q2;
    else if (q2 < 9.)
      rpigg = AEM_3PI * (16.32 + 2. * log_q2) + 0.00238 * log1p(3.927 * q2);
    else if (q2 < 1.e4)
      rpigg = AEM_3PI * (13.4955 + 3. * log_q2) + 0.00165 + 0.00299 * log_1_pl_q2;
    else
      rpigg = AEM_3PI * (13.4955 + 3. * log_q2) + 0.00221 + 0.00293 * log_1_pl_q2;
    return constants::ALPHA_EM / (1. - rpigg);
  }

private:
  static constexpr double AEM_3PI = constants::ALPHA_EM / 3. * M_1_PI;
  const double q2min_;
};
REGISTER_ALPHAEM_MODULE("burkhardt", AlphaEMBurkhardt);
