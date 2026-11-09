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

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGenHerwig6/DiffractiveSystemSeed.h"

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
  }  // namespace
}  // namespace cepgen::herwig6

using namespace cepgen::herwig6;

DiffractiveSystemSeed DiffractiveSystemSeed::build(const Momentum& beam,
                                                   const Momentum& system,
                                                   long beam_pdg_id,
                                                   const double* rmass,
                                                   utils::RandomGenerator& rng,
                                                   double primordial_kt) {
  if (std::labs(beam_pdg_id) != 2212)
    throw CG_FATAL("herwig6:DiffractiveSystemSeed")
        << "Beam particle with PDG id " << beam_pdg_id
        << " is not supported. Only protons and antiprotons can be dissociated.";

  DiffractiveSystemSeed diffractive_seed;
  diffractive_seed.system_mass = system.mass();

  // valence flavour content, weighted by the quark electric charge squared and its multiplicity in the proton
  const bool struck_down = rng.uniform(0., 1.) < 1. / 9.;
  const bool anti = beam_pdg_id < 0;
  diffractive_seed.struck_id =
      struck_down ? (anti ? InternalId::dBar : InternalId::down) : (anti ? InternalId::uBar : InternalId::up);
  diffractive_seed.spectator_id =
      struck_down ? (anti ? InternalId::uu_bar : InternalId::uu) : (anti ? InternalId::ud0_bar : InternalId::ud0);
  diffractive_seed.struck_mass = rmass[static_cast<int>(diffractive_seed.struck_id)];
  diffractive_seed.spectator_mass = rmass[static_cast<int>(diffractive_seed.spectator_id)];

  const auto mass2 = diffractive_seed.system_mass * diffractive_seed.system_mass;
  const auto sum = diffractive_seed.struck_mass + diffractive_seed.spectator_mass,
             diff = diffractive_seed.struck_mass - diffractive_seed.spectator_mass;
  const auto lambda = (mass2 - sum * sum) * (mass2 - diff * diff);
  if (lambda <= 0.)
    throw CG_FATAL("herwig6:DiffractiveSystemSeed")
        << "Dissociated system of mass " << diffractive_seed.system_mass << " GeV is too light for its constituents.";
  const auto p_star = 0.5 * std::sqrt(lambda) / diffractive_seed.system_mass;

  // the spectator follows the beam direction, as seen from the rest frame of the dissociated system
  const auto beam_in_system =
      Momentum(beam).lorentzBoost(Momentum::fromPxPyPzE(-system.px(), -system.py(), -system.pz(), system.energy()));
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
    diffractive_seed.spectator_momentum[i] = p_star * direction[i];
  return diffractive_seed;
}
