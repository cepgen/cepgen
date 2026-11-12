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

#include <CepGen/Core/Exception.h>
#include <CepGen/Event/Particle.h>
#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/PDG.h>
#include <CepGen/Utils/RandomGenerator.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "CepGenHerwig6/DiffractiveState.h"
#include "CepGenHerwig6/Herwig6Interface.h"

namespace cepgen::herwig6 {
  namespace {
    constexpr int kBeamStatus = 103;       // spectator entry holding the system kinematics
    constexpr int kPartonStatus = 149;     // final-state parton, copied to status 150 by HWCGSP
    constexpr int kSpectatorStatus = 147;  // beam remnant, copied to status 153 by HWCGSP

    void fillParton(
        int index, int status, DiffractiveState::InternalId id, const Momentum& p3, double mass, int partner) {
      hepevt_.isthep[index] = status;
      hepevt_.idhep[index] = hwprop_.idpdg[static_cast<int>(id)];
      hwevnt_.idhw[index] = static_cast<int>(id);
      hepevt_.jmohep[index][0] = 1;
      // a colour singlet pair is connected through both of its colour lines
      hepevt_.jmohep[index][1] = partner + 1;
      hepevt_.jdahep[index][1] = partner + 1;
      hepevt_.phep[index][0] = p3.px();
      hepevt_.phep[index][1] = p3.py();
      hepevt_.phep[index][2] = p3.pz();
      hepevt_.phep[index][3] = std::hypot(p3.p(), mass);
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
}  // namespace cepgen::herwig6

using namespace cepgen::herwig6;

DiffractiveState::DiffractiveState(utils::RandomGenerator& rng,
                                   const Particle& beam,
                                   const Momentum& system,
                                   double primordial_kt)
    : hadr::DiffractiveState(rng, beam.integerPdgId(), {{PDG::down, PDG::uu1, 1.}, {PDG::up, PDG::ud0, 8.}}),
      valquark_id_(internalId(valquark_pdg_id_)),
      diquark_id_(internalId(diquark_pdg_id_)) {
  initialise(system, beam.momentum(), primordial_kt);
}

bool DiffractiveState::hadronise(hadr::Fragments& fragments, double tolerance) const {
  fragments.clear();
  hwuine();

  // the system is hadronised at rest, as a post-shower colour singlet
  hepevt_.nhep = 3;
  hepevt_.isthep[0] = kBeamStatus;
  hepevt_.idhep[0] = 0;
  hwevnt_.idhw[0] = static_cast<int>(DiffractiveState::InternalId::centreOfMass);
  hepevt_.phep[0][0] = hepevt_.phep[0][1] = hepevt_.phep[0][2] = 0.;
  hepevt_.phep[0][3] = hepevt_.phep[0][4] = system_mass_;
  fillParton(1, kPartonStatus, valquark_id_, -diquark_momentum_, valquark_mass_, 2);
  fillParton(2, kSpectatorStatus, diquark_id_, diquark_momentum_, diquark_mass_, 1);

  hwcfor();
  hwcdec();
  hwdhad();
  hwdhvy();
  hwufne();  // mandatory: HWUINE aborts the run if it was not called for the previous event
  if (hwevnt_.ierror != 0)
    return false;

  for (int i = 0; i < hepevt_.nhep; ++i) {
    const auto status = hepevt_.isthep[i];
    const bool stable = status == 1;
    if (!stable && (status < 190 || status > 200))  // only keep hadrons
      continue;
    registerParticle(hwevnt_.idhw[i], hepevt_.idhep[i]);
    fragments.emplace_back(hadr::Fragment{
        i,
        hepevt_.jmohep[i][0] - 1,
        hepevt_.idhep[i],
        stable,
        Momentum::fromPxPyPzE(hepevt_.phep[i][0], hepevt_.phep[i][1], hepevt_.phep[i][2], hepevt_.phep[i][3])});
  }
  return momentumConserved(fragments, tolerance);
}

double DiffractiveState::massFromPdgId(long pdg_id) const {
  return hwprop_.rmass[static_cast<int>(internalId(pdg_id))];
}

DiffractiveState::InternalId DiffractiveState::internalId(long pdg_id) {
  const bool anti = pdg_id < 0;
  switch (std::labs(pdg_id)) {
    case PDG::down:
      return anti ? InternalId::down_bar : InternalId::down;
    case PDG::up:
      return anti ? InternalId::up_bar : InternalId::up;
    case PDG::uu1:
      return anti ? InternalId::uu_bar : InternalId::uu;
    case PDG::ud0:
      return anti ? InternalId::ud0_bar : InternalId::ud0;
    default:
      throw CG_FATAL("herwig6:DiffractiveState") << "No Herwig 6 identifier for parton with PDG id " << pdg_id << ".";
  }
}
