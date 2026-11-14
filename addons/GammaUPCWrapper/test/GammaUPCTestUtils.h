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

#ifndef CepGenGammaUPC_test_GammaUPCTestUtils_h
#define CepGenGammaUPC_test_GammaUPCTestUtils_h

#include <CepGen/Core/Exception.h>
#include <CepGen/Event/Event.h>
#include <CepGen/Physics/HeavyIon.h>
#include <CepGen/Physics/PDG.h>

#include <cmath>
#include <string>

#include "CepGenGammaUPC/GammaUPCInterface.h"

namespace gammaUPC::test {
  constexpr double kMassProton = 0.938272081;  ///< proton mass used internally by gammaUPC (GeV)
  constexpr double kMassNucleon = 0.9315;      ///< average nucleon mass in nuclei used internally by gammaUPC (GeV)

  /// Beam definition, with the full-ion kinematics convention (E = A * E_N, M = M_A) used in the event content
  struct BeamSpec {
    static BeamSpec proton(double energy) {
      return BeamSpec{cepgen::PDG::proton, energy, cepgen::PDG::get().mass(cepgen::PDG::proton), energy, 1, 1};
    }
    static BeamSpec ion(unsigned short a, cepgen::Element z, double energy_per_nucleon) {
      const cepgen::HeavyIon hi(a, z);
      return BeamSpec{static_cast<cepgen::pdgid_t>(hi),
                      a * energy_per_nucleon,
                      hi.mass(),
                      energy_per_nucleon,
                      a,
                      static_cast<int>(z)};
    }
    cepgen::pdgid_t pdg_id;
    double energy;              ///< total beam energy (GeV)
    double mass;                ///< total beam particle mass (GeV)
    double energy_per_nucleon;  ///< (GeV)
    int a, z;
    inline bool heavyIon() const { return a > 1; }
    /// Lorentz factor as computed internally by gammaUPC
    inline double gammaUPCLorentzFactor() const {
      return energy_per_nucleon / (heavyIon() ? kMassNucleon : kMassProton);
    }
  };

  /// Minimal elastic A B -> A (gamma gamma -> X) B event, with photon energies omega1 (emitted by the +z beam) and
  /// omega2 (emitted by the -z beam), carried by a two-particle central system (sum(E+pz) = 2 omega1,
  /// sum(E-pz) = 2 omega2); the photons per-nucleon momentum fractions are omega(1,2) / E_N(1,2)
  inline cepgen::Event buildEvent(const BeamSpec& beam1, const BeamSpec& beam2, double omega1, double omega2) {
    auto event = cepgen::Event::minimal(2);
    const auto pz1 = std::sqrt(beam1.energy * beam1.energy - beam1.mass * beam1.mass),
               pz2 = std::sqrt(beam2.energy * beam2.energy - beam2.mass * beam2.mass);
    event.oneWithRole(cepgen::Particle::Role::IncomingBeam1)
        .setIntegerPdgId(static_cast<long>(beam1.pdg_id))
        .setMomentum(0., 0., +pz1, beam1.energy);
    event.oneWithRole(cepgen::Particle::Role::IncomingBeam2)
        .setIntegerPdgId(static_cast<long>(beam2.pdg_id))
        .setMomentum(0., 0., -pz2, beam2.energy);
    auto central = event[cepgen::Particle::Role::CentralSystem];
    central.at(0).get().setMomentum(0., 0., +omega1, omega1);
    central.at(1).get().setMomentum(0., 0., -omega2, omega2);
    return event;
  }

  /// Set the gammaUPC beams configuration (total beam energies, (Z, A) content), positive-z beam first
  inline BeamsContent configure(const BeamSpec& beam1, const BeamSpec& beam2) {
    to_collider_.ebeam = {beam1.energy, beam2.energy};
    to_collider_.lpp = {2, 2};
    return setBeamsContent({beam1.z, beam2.z}, {beam1.a, beam2.a});
  }

  /// Within its scope, fatal exceptions do not terminate the process once caught (e.g. to test them)
  class ScopedNonTerminatingFatals {
  public:
    ScopedNonTerminatingFatals() : previous_handler_(cepgen::Exception::setFatalHandler([] {})) {}
    ~ScopedNonTerminatingFatals() { cepgen::Exception::setFatalHandler(previous_handler_); }

  private:
    const cepgen::Exception::FatalHandler previous_handler_;
  };

  inline std::string label(double v1, double v2) { return "(" + std::to_string(v1) + ", " + std::to_string(v2) + ")"; }
}  // namespace gammaUPC::test

#endif
