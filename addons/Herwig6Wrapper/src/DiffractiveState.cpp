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
    using Vec3 = std::array<double, 3>;
    Vec3 cross(const Vec3& a, const Vec3& b) {
      return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    }
    double norm(const Vec3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
    Vec3 unit(const Vec3& a) {
      const auto n = norm(a);
      return {a[0] / n, a[1] / n, a[2] / n};
    }
    Vec3 randomDirection(utils::RandomGenerator& rng) {
      const auto cos_theta = rng.uniform(-1., 1.), phi = rng.uniform(0., 2. * M_PI);
      const auto sin_theta = std::sqrt(1. - cos_theta * cos_theta);
      return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
    }
    constexpr int kBeamStatus = 103;       // spectator entry holding the system kinematics
    constexpr int kPartonStatus = 149;     // final-state parton, copied to status 150 by HWCGSP
    constexpr int kSpectatorStatus = 147;  // beam remnant, copied to status 153 by HWCGSP

    void fillParton(int index, int status, DiffractiveState::InternalId id, const double* p3, double mass, int partner) {
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
}  // namespace cepgen::herwig6

using namespace cepgen::herwig6;

DiffractiveState::DiffractiveState(const Particle& beam,
                                   const Momentum& system,
                                   const double* rmass,
                                   utils::RandomGenerator& rng,
                                   double primordial_kt)
    : hadr::Hadroniser::DiffractiveState(system) {
  const auto beam_pdg_id = beam.integerPdgId();
  if (std::labs(beam_pdg_id) != PDG::proton)
    throw CG_FATAL("herwig6:DiffractiveState") << "Beam particle with PDG id " << beam.pdgId()
                                               << " is not supported. Only protons and antiprotons can be dissociated.";

  // valence flavour content, weighted by the quark electric charge squared and its multiplicity in the proton
  const bool valquark_down = rng.uniform(0., 1.) < 1. / 9.;
  if (beam_pdg_id >= 0) {  // quark/diquark
    valquark_id_ = valquark_down ? InternalId::down : InternalId::up;
    diquark_id_ = valquark_down ? InternalId::uu : InternalId::ud0;
  } else {  // anti-quark/anti-diquark
    valquark_id_ = valquark_down ? InternalId::down_bar : InternalId::up_bar;
    diquark_id_ = valquark_down ? InternalId::uu_bar : InternalId::ud0_bar;
  }
  valquark_mass_ = rmass[static_cast<int>(valquark_id_)];
  diquark_mass_ = rmass[static_cast<int>(diquark_id_)];

  const auto mass2 = system_mass_ * system_mass_;
  const auto sum = valquark_mass_ + diquark_mass_, diff = valquark_mass_ - diquark_mass_;
  const auto lambda = (mass2 - sum * sum) * (mass2 - diff * diff);
  if (lambda <= 0.)
    throw CG_FATAL("herwig6:DiffractiveState")
        << "Dissociated system of mass " << system_mass_ << " GeV is too light for its constituents.";
  const auto p_star = 0.5 * std::sqrt(lambda) / system_mass_;

  // the spectator follows the beam direction, as seen from the rest frame of the dissociated system
  const auto beam_in_system =
      Momentum(beam.momentum())
          .lorentzBoost(Momentum::fromPxPyPzE(-system.px(), -system.py(), -system.pz(), system.energy()));
  const Vec3 beam_vector{beam_in_system.px(), beam_in_system.py(), beam_in_system.pz()};
  const auto axis = norm(beam_vector) > 0. ? unit(beam_vector) : randomDirection(rng);

  Vec3 direction = axis;
  if (primordial_kt > 0.) {
    const auto px = rng.normal(0., primordial_kt * M_SQRT1_2), py = rng.normal(0., primordial_kt * M_SQRT1_2);
    const auto pt = std::hypot(px, py);
    if (pt > 0.) {
      const auto sin_theta = std::min(pt / p_star, 1.), cos_theta = std::sqrt(1. - sin_theta * sin_theta);
      const Vec3 helper = std::fabs(axis[0]) < 0.9 ? Vec3{1., 0., 0.} : Vec3{0., 1., 0.};
      const auto e1 = unit(cross(axis, helper)), e2 = cross(axis, e1);
      for (size_t i = 0; i < 3; ++i)
        direction[i] = cos_theta * axis[i] + sin_theta * (px / pt * e1[i] + py / pt * e2[i]);
    }
  }
  for (size_t i = 0; i < 3; ++i)
    diquark_momentum_[i] = p_star * direction[i];
}

bool DiffractiveState::hadronise(hadr::Hadroniser::Fragments& fragments, double tolerance) const {
  fragments.clear();
  hwuine();

  // the system is hadronised at rest, as a post-shower colour singlet
  hepevt_.nhep = 3;
  hepevt_.isthep[0] = kBeamStatus;
  hepevt_.idhep[0] = 0;
  hwevnt_.idhw[0] = static_cast<int>(DiffractiveState::InternalId::centreOfMass);
  hepevt_.phep[0][0] = hepevt_.phep[0][1] = hepevt_.phep[0][2] = 0.;
  hepevt_.phep[0][3] = hepevt_.phep[0][4] = system_mass_;
  const auto& spectator = diquark_momentum_;
  const double struck[3] = {-spectator[0], -spectator[1], -spectator[2]};
  fillParton(1, kPartonStatus, valquark_id_, struck, valquark_mass_, 2);
  fillParton(2, kSpectatorStatus, diquark_id_, spectator.data(), diquark_mass_, 1);

  hwcfor();
  hwcdec();
  hwdhad();
  hwdhvy();
  hwufne();  // mandatory: HWUINE aborts the run if it was not called for the previous event
  if (hwevnt_.ierror != 0)
    return false;

  Momentum sum{0., 0., 0., system_mass_};
  for (int i = 0; i < hepevt_.nhep; ++i) {
    const auto status = hepevt_.isthep[i];
    const bool stable = status == 1;
    if (!stable && (status < 190 || status > 200))  // only keep hadrons
      continue;
    registerParticle(hwevnt_.idhw[i], hepevt_.idhep[i]);
    auto& fragment = fragments.emplace_back(hadr::Hadroniser::Fragment{
        i,
        hepevt_.jmohep[i][0] - 1,
        hepevt_.idhep[i],
        stable,
        Momentum::fromPxPyPzE(hepevt_.phep[i][0], hepevt_.phep[i][1], hepevt_.phep[i][2], hepevt_.phep[i][3])});
    if (fragment.stable)
      sum -= fragment.momentum;
  }
  if (const auto residual = sum.p(); residual > tolerance * system_mass_) {
    CG_WARNING("herwig6:hadronise") << "Momentum not conserved down to " << tolerance
                                    << " tolerance after hadronisation (residual=" << residual
                                    << " GeV for a system of mass " << system_mass_ << " GeV).";
    return false;
  }
  return true;
}
