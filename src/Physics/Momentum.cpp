/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2015-2025  Laurent Forthomme
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

#include <cmath>
#include <iomanip>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/Momentum.h"
#include "CepGen/Utils/Algebra.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/Math.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;

/// Express an angle in between two extrema
static double normalisePhi(double phi, const Limits& range) {
  static constexpr double M_2PI = 2. * M_PI;
  if (range.range() != M_2PI)
    throw CG_FATAL("Momentum:normalisePhi") << "Invalid boundaries for the angle normalisation: " << range
                                            << ". Must be of range 2*pi, has range " << range.range() << ".";
  while (phi < range.min())
    phi += M_2PI;
  while (phi > range.max())
    phi -= M_2PI;
  return phi;  // contained in range
}

static double normaliseSqrt(double x2) { return std::sqrt(x2 < 0. ? -x2 : x2); }

Momentum::Momentum(double x, double y, double z, double t) : std::array<double, 4>{{x, y, z, t == -1. ? 0. : t}} {
  computeP();
}

Momentum::Momentum(double* p) {
  std::copy_n(p, 4, begin());
  computeP();
}

Momentum::Momentum(const Vector& vec) {
  if (vec.size() < 3)
    throw CG_FATAL("Momentum") << "Failed to initialise a momentum from a vector with coordinates " << vec
                               << ". Should have at least 3 coordinates.";
  setPx(vec(0)).setPy(vec(1)).setPz(vec(2));
  if (vec.size() > 3)
    setEnergy(vec(3));
}

bool Momentum::operator==(const Momentum& oth) const {
  return px() == oth.px() && py() == oth.py() && pz() == oth.pz() && energy() == oth.energy();
}

//--- static constructors

Momentum Momentum::fromPtEtaPhiE(double pt, double eta, double phi, double e) {
  return Momentum(pt * cos(phi), pt * sin(phi), pt * sinh(eta), e);
}

Momentum Momentum::fromPtEtaPhiM(double pt, double eta, double phi, double m) {
  return Momentum::fromPtEtaPhiE(pt, eta, phi, utils::fastHypot(pt * cosh(eta), m));
}

Momentum Momentum::fromPThetaPhiE(double p, double theta, double phi, double e) {
  const double pt = p * sin(theta), px = pt * cos(phi), py = pt * sin(phi);
  return Momentum(px, py, p * cos(theta), e);
}

Momentum Momentum::fromPxPyPzE(double px, double py, double pz, double e) { return Momentum(px, py, pz, e); }

Momentum Momentum::fromPxPyPzM(double px, double py, double pz, double m) {
  return Momentum(px, py, pz).setMass(m).computeP();
}

Momentum Momentum::fromPxPyYM(double px, double py, double rap, double m) {
  const double et = utils::fastHypot(px, py, m);
  return Momentum(px, py, et * sinh(rap), et * cosh(rap));
}

Momentum Momentum::fromPtYPhiE(double pt, double rap, double phi, double e) {
  return Momentum(pt * cos(phi), pt * sin(phi), e * tanh(rap), e);
}

Momentum Momentum::fromPtYPhiM(double pt, double rap, double phi, double m) {
  const double et = utils::fastHypot(pt, m);
  return Momentum(pt * cos(phi), pt * sin(phi), et * sinh(rap), et * cosh(rap));
}

//--- arithmetic operators

// unary operators

Momentum Momentum::operator+() const { return Momentum(*this); }

Momentum Momentum::operator-() const { return Momentum(-px(), -py(), -pz(), energy()); }

// operators with scalars

Momentum Momentum::operator*(double c) const { return Momentum(c * px(), c * py(), c * pz(), c * energy()); }

Momentum& Momentum::operator*=(double c) {
  *this = *this * c;
  return computeP();
}

Momentum Momentum::operator/(double c) const { return Momentum(px() / c, py() / c, pz() / c, energy() / c); }

Momentum& Momentum::operator/=(double c) {
  *this = *this / c;
  return computeP();
}

// operators with other momenta

Momentum Momentum::operator+(const Momentum& mom) const {
  return Momentum(px() + mom.px(), py() + mom.py(), pz() + mom.pz(), energy() + mom.energy());
}

Momentum& Momentum::operator+=(const Momentum& mom) {
  *this = *this + mom;
  return computeP();
}

Momentum Momentum::operator-(const Momentum& mom) const {
  return Momentum(px() - mom.px(), py() - mom.py(), pz() - mom.pz(), energy() - mom.energy());
}

Momentum& Momentum::operator-=(const Momentum& mom) {
  *this = *this - mom;
  return computeP();
}

double Momentum::threeProduct(const Momentum& mom) const {
  CG_DEBUG_LOOP("Momentum") << "  (" << px() << ", " << py() << ", " << pz() << ")\n\t"
                            << "* (" << mom.px() << ", " << mom.py() << ", " << mom.pz() << ")\n\t"
                            << "= " << px() * mom.px() + py() * mom.py() + pz() * mom.pz();
  return px() * mom.px() + py() * mom.py() + pz() * mom.pz();
}

double Momentum::fourProduct(const Momentum& mom) const {
  CG_DEBUG_LOOP("Momentum") << "  (" << px() << ", " << py() << ", " << pz() << ", " << energy() << ")\n\t"
                            << "* (" << mom.px() << ", " << mom.py() << ", " << mom.pz() << ", " << mom.energy()
                            << ")\n\t"
                            << "= " << energy() * mom.energy() - threeProduct(mom);
  return energy() * mom.energy() - threeProduct(mom);
}

double Momentum::crossProduct(const Momentum& mom, coord_t direction) const {
  switch (direction) {
    case X:
      return py() * mom.pz() - pz() * mom.py();
    case Y:
      return px() * mom.pz() - pz() * mom.px();
    case Z:
      return px() * mom.py() - py() * mom.px();
    case E:
    default:
      throw CG_FATAL("Momentum:crossProduct")
          << "Unhandled coordinate for the cross-product: only (x, y, z) spatial coordinates are supported.";
  }
}

double Momentum::operator*(const Momentum& mom) const { return threeProduct(mom); }

Momentum Momentum::operator%(const Momentum& mom) const {
  return Momentum(crossProduct(mom, X), crossProduct(mom, Y), crossProduct(mom, Z));
}

//--- various setters

Momentum& Momentum::setPx(double px) {
  (*this)[X] = px;
  return computeP();
}

Momentum& Momentum::setPy(double py) {
  (*this)[Y] = py;
  return computeP();
}

Momentum& Momentum::setPz(double pz) {
  (*this)[Z] = pz;
  return computeP();
}

Momentum& Momentum::setEnergy(double e) {
  (*this)[E] = e;
  return *this;
}

Momentum& Momentum::setMass(double m) { return setEnergy(utils::fastHypot(p_, m)).computeP(); }

Momentum& Momentum::setMass2(double m2) { return setEnergy(std::sqrt(p_ * p_ + m2)).computeP(); }

Momentum& Momentum::setP(double px, double py, double pz, double e) { return setEnergy(e).setP(px, py, pz); }

Momentum& Momentum::setP(double px, double py, double pz) { return setPx(px).setPy(py).setPz(pz); }

//--- kinematics constraints

Momentum& Momentum::computeEnergyFromMass(double on_shell_mass) { return setMass(on_shell_mass); }

Momentum& Momentum::computePzFromMass(double on_shell_mass) {
  return setPz(normaliseSqrt(pz() * pz() + mass2() - on_shell_mass * on_shell_mass));
}

Momentum& Momentum::computeP() {
  p_ = utils::fastHypot(px(), py(), pz());
  return *this;
}

Momentum& Momentum::truncate(double tolerance) {
  std::replace_if(begin(), end(), [&tolerance](const auto& p) { return p <= tolerance; }, 0.);
  return computeP();
}

//--- various getters

std::array<double, 5> Momentum::pVector() const {
  std::array<double, 5> out;
  std::copy(begin(), end(), out.begin());
  out[4] = mass();
  return out;
}

Momentum::operator Vector() const { return Vector{px(), py(), pz(), energy()}; }

double Momentum::energyT2() const {
  const auto square_pt_value = pt2();
  return square_pt_value > 0. ? energy2() * square_pt_value / (square_pt_value + pz() * pz()) : 0.;
}

double Momentum::energyT() const { return normaliseSqrt(energyT2()); }

double Momentum::mass2() const { return energy2() - p2(); }

double Momentum::mass() const { return normaliseSqrt(mass2()); }

double Momentum::massT2() const { return energy2() - pz() * pz(); }

double Momentum::massT() const { return normaliseSqrt(massT2()); }

double Momentum::theta() const { return atan2(pt(), pz()); }

double Momentum::phi() const { return normalisePhi(atan2(py(), px()), {0., 2. * M_PI}); }

double Momentum::pt() const { return utils::fastHypot(px(), py()); }

double Momentum::pt2() const { return px() * px() + py() * py(); }

Momentum Momentum::transverse(coord_t direction) const {
  switch (direction) {
    case X:
      return Momentum::fromPxPyPzE(0., py(), pz(), utils::fastSqrtSqDiff(energy(), px(), utils::Normalise::yes));
    case Y:
      return Momentum::fromPxPyPzE(px(), 0., pz(), utils::fastSqrtSqDiff(energy(), py(), utils::Normalise::yes));
    case Z:
      return Momentum::fromPxPyPzE(px(), py(), 0., utils::fastSqrtSqDiff(energy(), pz(), utils::Normalise::yes));
    default:
      throw CG_FATAL("Momentum:transverse")
          << "Invalid direction of interest for the transverse component computation: " << static_cast<int>(direction)
          << ".";
  }
}

Momentum Momentum::longitudinal(coord_t direction) const {
  switch (direction) {
    case X:
      return Momentum::fromPxPyPzE(px(), 0., 0., utils::fastSqrtSqDiff(energy(), py(), pz(), utils::Normalise::yes));
    case Y:
      return Momentum::fromPxPyPzE(0., py(), 0., utils::fastSqrtSqDiff(energy(), px(), pz(), utils::Normalise::yes));
    case Z:
      return Momentum::fromPxPyPzE(0., 0., pz(), utils::fastSqrtSqDiff(energy(), px(), py(), utils::Normalise::yes));
    default:
      throw CG_FATAL("Momentum:longitudinal")
          << "Invalid longitudinal direction: " << static_cast<int>(direction) << ".";
  }
}

double Momentum::eta() const {
  const auto pt_value = pt();
  return (utils::positive(pt_value) ? std::log((p() + std::fabs(pz())) / pt_value) : 9999.) * utils::sign(pz());
}

double Momentum::rapidity() const {
  return utils::positive(energy()) ? std::log((energy() + pz()) / (energy() - pz())) * 0.5
                                   : std::numeric_limits<double>::infinity() * utils::sign(pz());
}

double Momentum::deltaEta(const Momentum& oth) const { return std::fabs(eta() - oth.eta()); }

double Momentum::deltaPhi(const Momentum& oth) const {
  return normalisePhi(phi() - oth.phi(), {-M_PI, M_PI});  // has to be contained in [-M_PI, M_PI]
}

double Momentum::deltaPt(const Momentum& oth) const { return std::fabs(pt() - oth.pt()); }

double Momentum::deltaR(const Momentum& oth) const {
  return utils::fastHypot(rapidity() - oth.rapidity(), deltaPhi(oth));
}

//--- boosts/rotations

double Momentum::beta() const {
  const auto mom = p(), ene = energy();
  if (ene == 0.) {
    if (mom == 0.)
      return 0.;
    CG_WARNING("Momentum:beta") << "beta computed for t=0 momentum.";
    return 1. / ene;
  }
  if (mass2() <= 0.)
    CG_WARNING("Momentum:beta") << "beta computed for an invalid, non-time-like momentum.";
  return mom / ene;
}

double Momentum::gamma2() const {
  const auto squared_momentum = p2(), squared_energy = energy2();
  if (squared_energy == 0.) {
    if (squared_momentum == 0.)
      return 1.;
    CG_WARNING("Momentum:gamma") << "gamma computed for t=0 momentum.";
  }
  if (squared_energy < squared_momentum) {
    CG_WARNING("Momentum:gamma") << "gamma computed for an invalid space-like momentum.";
    return 0.;
  }
  if (squared_energy == squared_momentum)
    CG_WARNING("Momentum:gamma") << "gamma computed for a light-like momentum.";
  return squared_energy / (squared_energy - squared_momentum);
}

double Momentum::gamma() const { return std::sqrt(gamma2()); }

Momentum& Momentum::betaGammaBoost(double gamma, double beta_gamma) {
  if (gamma == 1. && beta_gamma == 0.)
    return *this;  // trivial case
  const auto apz = pz(), ae = energy();
  return setEnergy(gamma * ae + beta_gamma * apz).setPz(gamma * apz + beta_gamma * ae);
}

Momentum& Momentum::lorentzBoost(const Momentum& mom) {
  if (mom.p() == 0.)
    return *this;  // do not boost on a system at rest
  const auto mass = mom.mass();
  if (mass == 0.)
    return *this;
  const auto pf4 = (threeProduct(mom) + energy() * mom.energy()) / mass;
  const auto fn = (pf4 + energy()) / (mass + mom.energy());
  (*this) += fn * mom;
  return setEnergy(pf4);
  /*const auto norm_mom = mom * (1. / mom.energy());
    const auto b2 = norm_mom.p2();
    const auto gamma = std::sqrt(1. / (1. - b2)), gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.;
    const auto bp = threeProduct(norm_mom);
    return setPx(px() + gamma2 * bp * norm_mom.px() + gamma * norm_mom.px() * energy())
        .setPy(py() + gamma2 * bp * norm_mom.py() + gamma * norm_mom.py() * energy())
        .setPz(pz() + gamma2 * bp * norm_mom.pz() + gamma * norm_mom.pz() * energy())
        .setEnergy(gamma * (energy() + bp));*/
}

Momentum& Momentum::rotatePhi(double phi, double sign) {
  const auto sin_phi = std::sin(phi), cos_phi = std::cos(phi);
  const auto pxp = px() * cos_phi + sign * py() * sin_phi, pyp = -px() * sin_phi + sign * py() * cos_phi;
  return setPx(pxp).setPy(pyp);
}

Momentum& Momentum::rotateThetaPhi(double theta, double phi) {
  const double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
  const double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
  double rotation_matrix[3][3];  //FIXME check this! cos(phi)->-sin(phi) & sin(phi)->cos(phi) --> phi->phi+pi/2 ?
  auto mom = std::vector(3, 0.);
  rotation_matrix[X][X] = -sin_phi;
  rotation_matrix[X][Y] = -cos_theta * cos_phi;
  rotation_matrix[X][Z] = sin_theta * cos_phi;
  rotation_matrix[Y][X] = cos_phi;
  rotation_matrix[Y][Y] = -cos_theta * sin_phi;
  rotation_matrix[Y][Z] = sin_theta * sin_phi;
  rotation_matrix[Z][X] = 0.;
  rotation_matrix[Z][Y] = sin_theta;
  rotation_matrix[Z][Z] = cos_theta;
  for (size_t i = X; i <= Z; ++i)
    for (size_t j = X; j <= Z; ++j)
      mom[i] += rotation_matrix[i][j] * (*this)[j];
  return setP(mom.at(X), mom.at(Y), mom.at(Z));
}

namespace cepgen {
  Momentum operator*(double scalar, const Momentum& mom) {
    return Momentum(scalar * mom.px(), scalar * mom.py(), scalar * mom.pz(), scalar * mom.energy());
  }

  std::ostream& operator<<(std::ostream& os, const Momentum& mom) {
    return os << utils::format(
               "(%8.3f,%8.3f,%8.3f;%8.3f|%8.3f)", mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass());
  }
}  // namespace cepgen
