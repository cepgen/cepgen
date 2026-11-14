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
#include <CepGen/Generator.h>
#include <CepGen/Modules/PartonFluxFactory.h>
#include <CepGen/Modules/SurvivalFactorFactory.h>
#include <CepGen/PartonFluxes/CollinearFlux.h>
#include <CepGen/Physics/Beam.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/SurvivalFactors/SurvivalFactor.h>
#include <CepGen/Utils/ArgumentsParser.h>
#include <CepGen/Utils/Limits.h>
#include <CepGen/Utils/Test.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "CepGenGammaUPC/GammaUPCInterface.h"

using namespace cepgen;
using namespace std::string_literals;

namespace {
  using XPoints = std::vector<std::pair<double, double> >;
  std::string xLabel(double x1, double x2) { return "(" + std::to_string(x1) + ", " + std::to_string(x2) + ")"; }

  /// Minimal pp -> p(gamma gamma -> X)p event, with a central system carrying sum(E +/- pz) = x(1,2) sqrt(s)
  Event buildEvent(double beam_energy, double x1, double x2) {
    auto event = Event::minimal(2);
    const auto mp = PDG::get().mass(PDG::proton), pz = std::sqrt(beam_energy * beam_energy - mp * mp);
    event.oneWithRole(Particle::Role::IncomingBeam1).setIntegerPdgId(PDG::proton).setMomentum(0., 0., +pz, beam_energy);
    event.oneWithRole(Particle::Role::IncomingBeam2).setIntegerPdgId(PDG::proton).setMomentum(0., 0., -pz, beam_energy);
    const auto sqrt_s = event.cmEnergy();
    auto central = event[Particle::Role::CentralSystem];
    central.at(0).get().setMomentum(0., 0., +0.5 * x1 * sqrt_s, 0.5 * x1 * sqrt_s);  // E + pz = x1 sqrt(s)
    central.at(1).get().setMomentum(0., 0., -0.5 * x2 * sqrt_s, 0.5 * x2 * sqrt_s);  // E - pz = x2 sqrt(s)
    return event;
  }
}  // namespace

/// Unit tests for the gammaUPC wrappers (pp mode):
/// - bindings: the P_nohad=1 two-photon flux factorises into the single-photon (EDFF, no-b) fluxes;
/// - survival: the hadronic non-overlap probability reaches the Fortran code (S^2 strictly below 1, in a physical
///   range, decreasing with the two-photon mass; 1-5% effect on the pp luminosity, arXiv:2207.03012);
/// - "gammaUPC" collinear flux: once integrated over its ln(Q^2) range (as done by the collinear parton sampler),
///   it gives back x f(x), and vanishes outside;
/// - "gammaUPC" survival factor: retrieves x1/x2 from the central system, and reproduces the direct flux ratio.
int main(int argc, char* argv[]) {
  CG_TEST_EXIT_GUARD();  // gammaUPC may exit through Fortran STOP statements
  double sqrt_s, s2_min;
  std::vector<double> masses;
  bool debug;
  ArgumentsParser(argc, argv)
      .addOptionalArgument("sqrts,s", "two-proton centre-of-mass energy (GeV)", &sqrt_s, 14.e3)
      .addOptionalArgument(
          "masses,m", "two-photon masses to probe (GeV)", &masses, std::vector<double>{20., 100., 300.})
      .addOptionalArgument("s2min", "minimal acceptable survival factor", &s2_min, 0.85)
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .parse();
  CG_TEST_DEBUG(debug);
  initialise();

  // gammaUPC beams configuration (pp), identical to the one set by the survival factor module from the events, as
  // gammaUPC caches it (and its interpolation grids) at its first call
  const auto beam_energy = 0.5 * sqrt_s, gamma = beam_energy / PDG::get().mass(PDG::proton);
  to_collider_.ebeam = {beam_energy, beam_energy};
  to_collider_.lpp = {2, 2};
  // NOTE: never call it within a logging macro, as its stream may not be evaluated
  const auto beams_content = gammaUPC::setBeamsContent({1, 1}, {1, 1});
  CG_DEBUG("gammaupc_fluxes") << "Beams content passed through: " << beams_content << ".";

  const auto single_fluxes = [&gamma](double x1, double x2) {
    return gammaUPC::photonFluxP(x1, gamma) * gammaUPC::photonFluxP(x2, gamma);
  };
  const auto survival = [](double x1, double x2) {
    const auto flux_no_had = gammaUPC::twoPhotonFluxPP(x1, x2, 1);
    return flux_no_had > 0. ? gammaUPC::twoPhotonFluxPP(x1, x2, 0) / flux_no_had : -1.;
  };

  {  // bindings: P_nohad=1 two-photon flux = product of the single-photon fluxes
    CG_TEST_SET_PRECISION(1.e-5);
    for (const auto& [x1, x2] : XPoints{{1.e-3, 1.e-3}, {1.e-3, 2.e-2}, {5.e-2, 1.e-4}}) {
      const auto flux_no_had = gammaUPC::twoPhotonFluxPP(x1, x2, 1), product = single_fluxes(x1, x2);
      CG_TEST(flux_no_had > 0. && product > 0., "positive fluxes at x=" + xLabel(x1, x2));
      CG_TEST_EQUIV(flux_no_had / product, 1., "no-overlap two-photon flux factorisation at x=" + xLabel(x1, x2));
    }
    CG_TEST_RESET_PRECISION();
  }

  {  // survival probability for a central (y=0) two-photon system of mass m: x1 = x2 = m / sqrt(s)
    std::vector<double> s2_values;
    for (const auto& mass : masses) {
      const auto x = mass / sqrt_s, s2 = survival(x, x);
      CG_DEBUG("gammaupc_fluxes") << "m = " << mass << " GeV, x = " << x << ": S^2 = " << s2 << ".";
      CG_TEST(s2 >= s2_min && s2 < 1. - 1.e-6,
              "S^2 in [" + std::to_string(s2_min) + ", 1) at m=" + std::to_string(mass) +
                  " GeV (S^2=" + std::to_string(s2) + ")");
      s2_values.emplace_back(s2);
    }
    for (size_t i = 1; i < s2_values.size(); ++i)
      CG_TEST(s2_values.at(i) < s2_values.at(i - 1) + 1.e-3,
              "S^2 decreasing with mass (" + std::to_string(masses.at(i - 1)) + " -> " + std::to_string(masses.at(i)) +
                  " GeV)");
    CG_TEST_EQUIV(survival(2.e-3, 2.e-2), survival(2.e-2, 2.e-3), "pp S^2 symmetric under x1 <-> x2");
  }

  {  // collinear flux module
    const Limits q2_range{1.e-10, 5.};
    const auto log_q2_range = std::log(q2_range.max() / q2_range.min());
    auto flux = CollinearFluxFactory::get().build("gammaUPC", ParametersList().set("q2range", q2_range));
    CG_TEST(flux != nullptr, "gammaUPC collinear flux built");
    {  // fatal exceptions terminate the process once caught (by design): temporarily disable this behaviour
      const auto previous_handler = Exception::setFatalHandler([] {});
      CG_TEST_EXCEPT([&flux] { flux->fluxQ2(1.e-2, 1.); },
                     "collinear flux refuses to be evaluated before initialisation");
      Exception::setFatalHandler(previous_handler);
    }
    flux->initialise(
        Beam(ParametersList()
                 .set<int>("pdgId", PDG::proton)
                 .set<double>("pz", std::sqrt(beam_energy * beam_energy - std::pow(PDG::get().mass(PDG::proton), 2)))));
    CG_TEST_EQUAL(flux->partonPdgId(), (spdgid_t)PDG::photon, "collinear flux parton is a photon");
    CG_TEST_SET_PRECISION(1.e-6);
    for (const auto& x : {1.e-4, 1.e-3, 1.e-2, 1.e-1}) {
      const auto xf = x * gammaUPC::photonFluxP(x, gamma);
      // the collinear sampler integrates fluxQ2(x, t) / t dt = fluxQ2(x, t) d(ln t) over the q2 range
      double integral = 0.;
      const size_t num_steps = 100;
      for (size_t i = 0; i < num_steps; ++i) {  // midpoint rule in ln(Q^2)
        const auto log_q2 = std::log(q2_range.min()) + (i + 0.5) * log_q2_range / num_steps;
        integral += flux->fluxQ2(x, std::exp(log_q2));
      }
      integral *= log_q2_range / num_steps;
      CG_TEST_EQUIV(integral / xf, 1., "ln(Q^2)-integrated collinear flux = x f(x) at x=" + std::to_string(x));
    }
    CG_TEST_RESET_PRECISION();
    CG_TEST_EQUAL(flux->fluxQ2(1.e-2, 2. * q2_range.max()), 0., "collinear flux vanishes above the Q^2 range");
    CG_TEST_EQUAL(flux->fluxQ2(1.e-2, 0.5 * q2_range.min()), 0., "collinear flux vanishes below the Q^2 range");
  }

  {  // survival factor module
    auto survival_module =
        SurvivalFactorFactory::get().build("gammaUPC", ParametersList().set<bool>("forcePNonHadronisation", false));
    auto no_overlap_module =
        SurvivalFactorFactory::get().build("gammaUPC", ParametersList().set<bool>("forcePNonHadronisation", true));
    CG_TEST(survival_module && no_overlap_module, "gammaUPC survival factor modules built");
    CG_TEST_SET_PRECISION(1.e-5);
    for (const auto& [x1, x2] : XPoints{{1.e-3, 1.e-3}, {2.e-3, 2.e-2}, {1.e-2, 5.e-4}}) {
      const auto event = buildEvent(beam_energy, x1, x2);
      const auto s2 = (*survival_module)(event);
      CG_TEST_EQUIV(s2, survival(x1, x2), "module S^2 = direct flux ratio at x=" + xLabel(x1, x2));
      CG_TEST(s2 < 1., "module S^2 < 1 at x=" + xLabel(x1, x2));
      CG_TEST_EQUIV((*no_overlap_module)(event), 1., "module S^2 = 1 when forcing P_nohad=1, at x=" + xLabel(x1, x2));
    }
    CG_TEST_RESET_PRECISION();
    CG_TEST_EQUAL((*survival_module)(buildEvent(beam_energy, 0., 1.e-3)), 0., "module S^2 vanishes for x1 = 0");
  }

  CG_TEST_SUMMARY;
}
