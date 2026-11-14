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
#include <CepGen/Generator.h>
#include <CepGen/Modules/SurvivalFactorFactory.h>
#include <CepGen/SurvivalFactors/SurvivalFactor.h>
#include <CepGen/Utils/ArgumentsParser.h>
#include <CepGen/Utils/Test.h>

#include <vector>

#include "GammaUPCTestUtils.h"

using namespace cepgen;
using namespace gammaUPC::test;

/// Unit tests for the gammaUPC wrappers in ion-ion mode (Pb+Pb):
/// - bindings: the P_nohad=1 two-photon flux factorises into the nuclei single-photon fluxes, evaluated at the
///   per-nucleon momentum fractions and Lorentz factors;
/// - survival: S^2 = F / F(P_nohad=1) strictly below 1, in a physical range, decreasing with the two-photon mass,
///   symmetric under x1 <-> x2;
/// - "gammaUPC" survival factor: from an event with known photon energies (full-ion beam kinematics), reproduces the
///   direct S^2 at the per-nucleon x(1,2) = omega(1,2) / E_N
int main(int argc, char* argv[]) {
  CG_TEST_EXIT_GUARD();  // gammaUPC may exit through Fortran STOP statements
  double energy_n, s2_min, module_tolerance;
  std::vector<double> masses;
  bool debug;
  ArgumentsParser(argc, argv)
      .addOptionalArgument("energy-n", "lead beams energy per nucleon (GeV)", &energy_n, 2562.2)
      .addOptionalArgument("masses,m", "two-photon masses to probe (GeV)", &masses, std::vector<double>{10., 30., 100.})
      .addOptionalArgument("s2min", "minimal acceptable survival factor", &s2_min, 0.2)
      .addOptionalArgument("module-tolerance", "relative tolerance on module/direct S^2", &module_tolerance, 1.e-3)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();
  CG_TEST_DEBUG(debug);
  initialise();

  const auto lead = BeamSpec::ion(208, Element::Pb, energy_n);
  const auto mode = gammaUPC::HeavyIonMode::HardSphere;
  // NOTE: never call it within a logging macro, as its stream may not be evaluated
  const auto beams_content = configure(lead, lead);
  CG_DEBUG("gammaupc_fluxes_aa") << "Beams content passed through: " << beams_content << ".";
  CG_TEST(beams_content != gammaUPC::BeamsContent::none, "ions content passed to gammaUPC");

  const auto survival = [&mode](double x1, double x2) {
    const auto flux_no_had = gammaUPC::twoPhotonFluxAB(x1, x2, 1, mode);
    return flux_no_had > 0. ? gammaUPC::twoPhotonFluxAB(x1, x2, 0, mode) / flux_no_had : -1.;
  };

  {  // bindings: P_nohad=1 two-photon flux = nucleus flux x nucleus flux
    CG_TEST_SET_PRECISION(1.e-4);
    for (const auto& [x1, x2] :
         std::vector<std::pair<double, double> >{{1.e-3, 1.e-3}, {1.e-4, 1.e-2}, {2.e-2, 2.e-4}}) {
      const auto flux_no_had = gammaUPC::twoPhotonFluxAB(x1, x2, 1, mode),
                 product = gammaUPC::photonFluxA(x1, lead.gammaUPCLorentzFactor(), lead.z, 6.624 /* fm */) *
                           gammaUPC::photonFluxA(x2, lead.gammaUPCLorentzFactor(), lead.z, 6.624 /* fm */);
      CG_TEST(flux_no_had > 0. && product > 0., "positive fluxes at x=" + label(x1, x2));
      CG_TEST_EQUIV(flux_no_had / product, 1., "Pb-Pb no-overlap flux factorisation at x=" + label(x1, x2));
    }
    CG_TEST_RESET_PRECISION();
  }

  {  // survival probability for a central system of mass m: x1 = x2 = m / (2 E_N)
    std::vector<double> s2_values;
    for (const auto& mass : masses) {
      const auto x = 0.5 * mass / lead.energy_per_nucleon, s2 = survival(x, x);
      CG_DEBUG("gammaupc_fluxes_aa") << "m = " << mass << " GeV, x = " << x << ": S^2 = " << s2 << ".";
      CG_TEST(s2 >= s2_min && s2 < 1. - 1.e-6,
              "S^2 in [" + std::to_string(s2_min) + ", 1) at m=" + std::to_string(mass) +
                  " GeV (S^2=" + std::to_string(s2) + ")");
      s2_values.emplace_back(s2);
    }
    for (size_t i = 1; i < s2_values.size(); ++i)
      CG_TEST(s2_values.at(i) < s2_values.at(i - 1) + 1.e-3,
              "S^2 decreasing with mass (" + std::to_string(masses.at(i - 1)) + " -> " + std::to_string(masses.at(i)) +
                  " GeV)");
    CG_TEST_EQUIV(survival(1.e-3, 1.e-2), survival(1.e-2, 1.e-3), "Pb-Pb S^2 symmetric under x1 <-> x2");
  }

  {  // survival factor module
    auto module = SurvivalFactorFactory::get().build("gammaUPC",
                                                     ParametersList()
                                                         .set<bool>("forcePNonHadronisation", false)
                                                         .setAs<int, gammaUPC::HeavyIonMode>("heavyIonMode", mode));
    CG_TEST(module != nullptr, "gammaUPC survival factor module built");
    CG_TEST_SET_PRECISION(module_tolerance);
    for (const auto& [omega1, omega2] : std::vector<std::pair<double, double> >{{5., 5.}, {2., 40.}, {40., 2.}}) {
      const auto expected = survival(omega1 / lead.energy_per_nucleon, omega2 / lead.energy_per_nucleon),
                 s2 = (*module)(buildEvent(lead, lead, omega1, omega2));
      CG_DEBUG("gammaupc_fluxes_aa") << "omega = " << label(omega1, omega2) << " GeV: S^2(direct) = " << expected
                                     << ", S^2(module) = " << s2 << ".";
      CG_TEST_EQUIV(s2 / expected, 1., "Pb+Pb module S^2 = direct S^2 at omega=" + label(omega1, omega2) + " GeV");
    }
    CG_TEST_RESET_PRECISION();
  }

  CG_TEST_SUMMARY;
}
