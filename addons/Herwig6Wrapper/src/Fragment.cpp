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

#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/Message.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

#include "CepGenHerwig6/DiffractiveSystemSeed.h"
#include "CepGenHerwig6/Fragment.h"
#include "CepGenHerwig6/Herwig6Interface.h"

namespace cepgen::herwig6 {
  namespace {
    constexpr int kBeamStatus = 103;       // spectator entry holding the system kinematics
    constexpr int kPartonStatus = 149;     // final-state parton, copied to status 150 by HWCGSP
    constexpr int kSpectatorStatus = 147;  // beam remnant, copied to status 153 by HWCGSP
    constexpr double kMomentumTolerance = 1.e-6;

    void fillParton(
        int index, int status, DiffractiveSystemSeed::InternalId id, const double* p3, double mass, int partner) {
      hepevt_.isthep[index] = status;
      hepevt_.idhep[index] = hwprop_.idpdg[static_cast<int>(id)];
      hwevnt_.idhw[index] = static_cast<int>(id);
      hepevt_.jmohep[index][0] = 1;
      // a colour singlet pair is connected through both of its colour lines
      hepevt_.jmohep[index][1] = partner + 1;
      hepevt_.jdahep[index][1] = partner + 1;
      for (int i = 0; i < 3; ++i)
        hepevt_.phep[index][i] = p3[i];
      hepevt_.phep[index][3] = std::sqrt(p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2] + mass * mass);
      hepevt_.phep[index][4] = mass;
    }

    /// Make a Herwig 6 hadron known to the CepGen particles table
    void registerParticle(int id, long pdg_id) {
      if (PDG::get().has(std::labs(pdg_id)))
        return;
      std::string name(hwunam_.rname[id], 8);
      name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
      ParticleProperties properties;
      properties.pdgid = std::labs(pdg_id);
      properties.name = name;
      properties.human_name = name;
      properties.colours = 1;
      properties.mass = hwprop_.rmass[id];
      properties.width = hwprop_.rltim[id] > 0. ? constants::HBAR / hwprop_.rltim[id] : 0.;
      if (const auto charge = (pdg_id < 0 ? -1 : +1) * hwprop_.ichrg[id]; charge != 0)
        properties.charges = {charge, -charge};
      properties.fermion = std::fmod(hwprop_.rspin[id], 1.) != 0.;
      PDG::get().define(properties);
    }
  }  // namespace

  bool hadronise(const DiffractiveSystemSeed& seed, Fragments& fragments) {
    fragments.clear();
    hwuine();

    // the system is hadronised at rest, as a post-shower colour singlet
    const auto mass = seed.system_mass;
    hepevt_.nhep = 3;
    hepevt_.isthep[0] = kBeamStatus;
    hepevt_.idhep[0] = 0;
    hwevnt_.idhw[0] = static_cast<int>(DiffractiveSystemSeed::InternalId::centreOfMass);
    hepevt_.phep[0][0] = hepevt_.phep[0][1] = hepevt_.phep[0][2] = 0.;
    hepevt_.phep[0][3] = hepevt_.phep[0][4] = mass;
    const auto& spectator = seed.spectator_momentum;
    const double struck[3] = {-spectator[0], -spectator[1], -spectator[2]};
    fillParton(1, kPartonStatus, seed.struck_id, struck, seed.struck_mass, 2);
    fillParton(2, kSpectatorStatus, seed.spectator_id, spectator.data(), seed.spectator_mass, 1);

    hwcfor();
    hwcdec();
    hwdhad();
    hwdhvy();
    hwufne();  // mandatory: HWUINE aborts the run if it was not called for the previous event
    if (hwevnt_.ierror != 0)
      return false;

    double sum[4] = {0., 0., 0., 0.};
    for (int i = 0; i < hepevt_.nhep; ++i) {
      const auto status = hepevt_.isthep[i];
      const bool stable = status == 1;
      if (!stable && (status < 190 || status > 200))  // only keep hadrons
        continue;
      registerParticle(hwevnt_.idhw[i], hepevt_.idhep[i]);
      Fragment fragment;
      fragment.index = i;
      fragment.mother = hepevt_.jmohep[i][0] - 1;
      fragment.pdg_id = hepevt_.idhep[i];
      fragment.stable = stable;
      const auto* p = hepevt_.phep[i];
      fragment.momentum = Momentum::fromPxPyPzE(p[0], p[1], p[2], p[3]);
      fragments.emplace_back(fragment);
      if (stable)
        for (int j = 0; j < 4; ++j)
          sum[j] += p[j];
    }
    sum[3] -= mass;
    const auto residual = std::max({std::fabs(sum[0]), std::fabs(sum[1]), std::fabs(sum[2]), std::fabs(sum[3])});
    if (residual > kMomentumTolerance * mass) {
      CG_WARNING("herwig6:hadronise") << "Momentum not conserved down to " << kMomentumTolerance
                                      << " tolerance after hadronisation (residual=" << residual
                                      << " GeV for a system of mass " << mass << " GeV).";
      return false;
    }
    return true;
  }
}  // namespace cepgen::herwig6
