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

/// Unit tests for the gammaUPC wrappers in proton-ion mode (p+Pb, and its mirrored Pb+p configuration):
/// - bindings: the P_nohad=1 two-photon flux factorises into the proton and nucleus single-photon fluxes, evaluated at
///   the per-nucleon momentum fractions and Lorentz factors;
/// - survival: S^2 = F / F(P_nohad=1) strictly below 1, in a physical range, decreasing with the two-photon mass;
/// - "gammaUPC" survival factor: from an event with known photon energies (full-ion beam kinematics), reproduces the
///   direct S^2 at the per-nucleon x1 = omega1 / E_p, x2 = omega2 / E_N; and is invariant under the p+Pb -> Pb+p
///   mirroring of the event
int main(int argc, char* argv[]) {
  CG_TEST_EXIT_GUARD();  // gammaUPC may exit through Fortran STOP statements
  double energy_p, energy_n, s2_min, module_tolerance;
  std::vector<double> masses;
  bool debug;
  ArgumentsParser(argc, argv)
      .addOptionalArgument("energy-p", "proton beam energy (GeV)", &energy_p, 6500.)
      .addOptionalArgument("energy-n", "lead beam energy per nucleon (GeV)", &energy_n, 2562.2)
      .addOptionalArgument("masses,m", "two-photon masses to probe (GeV)", &masses, std::vector<double>{10., 50., 200.})
      .addOptionalArgument("s2min", "minimal acceptable survival factor", &s2_min, 0.5)
      .addOptionalArgument("module-tolerance", "relative tolerance on module/direct S^2", &module_tolerance, 1.e-3)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();
  CG_TEST_DEBUG(debug);
  initialise();

  const auto proton = BeamSpec::proton(energy_p), lead = BeamSpec::ion(208, Element::Pb, energy_n);
  const auto mode = gammaUPC::HeavyIonMode::HardSphere;
  // proton first: the configuration the survival factor module also sets for Pb+p (gammaUPC caches it at first call)
  CG_DEBUG("gammaupc_fluxes_pa") << "Beams content passed through: " << configure(proton, lead) << ".";

  const auto survival = [&mode](double x1, double x2) {  // x1: proton, x2: lead (per nucleon)
    const auto flux_no_had = gammaUPC::twoPhotonFluxPA(x1, x2, 1, mode);
    return flux_no_had > 0. ? gammaUPC::twoPhotonFluxPA(x1, x2, 0, mode) / flux_no_had : -1.;
  };

  {  // bindings: P_nohad=1 two-photon flux = proton flux x nucleus flux
    CG_TEST_SET_PRECISION(1.e-4);
    for (const auto& [x1, x2] :
         std::vector<std::pair<double, double> >{{1.e-3, 1.e-3}, {1.e-4, 1.e-2}, {2.e-2, 2.e-4}}) {
      const auto flux_no_had = gammaUPC::twoPhotonFluxPA(x1, x2, 1, mode),
                 product = gammaUPC::photonFluxP(x1, proton.gammaUPCLorentzFactor()) *
                           gammaUPC::photonFluxA(x2, lead.gammaUPCLorentzFactor(), lead.z, 6.624 /* fm */);
      CG_TEST(flux_no_had > 0. && product > 0., "positive fluxes at x=" + label(x1, x2));
      CG_TEST_EQUIV(flux_no_had / product, 1., "p-Pb no-overlap flux factorisation at x=" + label(x1, x2));
    }
    CG_TEST_RESET_PRECISION();
  }

  {  // survival probability for omega1 = omega2 = m / 2
    std::vector<double> s2_values;
    for (const auto& mass : masses) {
      const auto x1 = 0.5 * mass / proton.energy_per_nucleon, x2 = 0.5 * mass / lead.energy_per_nucleon,
                 s2 = survival(x1, x2);
      CG_DEBUG("gammaupc_fluxes_pa") << "m = " << mass << " GeV, x = " << label(x1, x2) << ": S^2 = " << s2 << ".";
      CG_TEST(s2 >= s2_min && s2 < 1. - 1.e-6,
              "S^2 in [" + std::to_string(s2_min) + ", 1) at m=" + std::to_string(mass) +
                  " GeV (S^2=" + std::to_string(s2) + ")");
      s2_values.emplace_back(s2);
    }
    for (size_t i = 1; i < s2_values.size(); ++i)
      CG_TEST(s2_values.at(i) < s2_values.at(i - 1) + 1.e-3,
              "S^2 decreasing with mass (" + std::to_string(masses.at(i - 1)) + " -> " + std::to_string(masses.at(i)) +
                  " GeV)");
  }

  {  // survival factor module
    const auto module_params = ParametersList()
                                   .set<bool>("forcePNonHadronisation", false)
                                   .setAs<int, gammaUPC::HeavyIonMode>("heavyIonMode", mode);
    auto module_pa = SurvivalFactorFactory::get().build("gammaUPC", module_params),
         module_ap = SurvivalFactorFactory::get().build("gammaUPC", module_params);  // each initialised at 1st event
    CG_TEST(module_pa && module_ap, "gammaUPC survival factor modules built");
    CG_TEST_SET_PRECISION(module_tolerance);
    for (const auto& [omega1, omega2] : std::vector<std::pair<double, double> >{{5., 5.}, {2., 40.}, {60., 3.}}) {
      const auto expected = survival(omega1 / proton.energy_per_nucleon, omega2 / lead.energy_per_nucleon);
      const auto s2_pa = (*module_pa)(buildEvent(proton, lead, omega1, omega2)),
                 s2_ap = (*module_ap)(buildEvent(lead, proton, omega2, omega1));  // mirrored event
      CG_DEBUG("gammaupc_fluxes_pa") << "omega = " << label(omega1, omega2) << " GeV: S^2(direct) = " << expected
                                     << ", S^2(p+Pb) = " << s2_pa << ", S^2(Pb+p) = " << s2_ap << ".";
      CG_TEST_EQUIV(s2_pa / expected, 1., "p+Pb module S^2 = direct S^2 at omega=" + label(omega1, omega2) + " GeV");
      CG_TEST_EQUIV(s2_ap / expected, 1., "Pb+p module S^2 = direct S^2 at omega=" + label(omega1, omega2) + " GeV");
    }
    CG_TEST_RESET_PRECISION();
  }

  CG_TEST_SUMMARY;
}
