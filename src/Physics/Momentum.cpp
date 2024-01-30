/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

namespace cepgen {
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

  const Matrix metricMatrix{
      {-1., 0., 0., 0.},  // x
      {0., -1., 0., 0.},  // y
      {0., 0., -1., 0.},  // z
      {0., 0., 0., 1.}    // t
  };

  Momentum::Momentum(double x, double y, double z, double t) : Vector{x, y, z, t == -1. ? 0. : t} { computeP(); }

  Momentum::Momentum(double* p) : Vector{p[0], p[1], p[2], p[3]} { computeP(); }

  Momentum::Momentum(const Matrix& mat) : Vector(mat.column(0)) { computeP(); }

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

  double Momentum::threeProduct(const Momentum& mom) const {
    CG_DEBUG_LOOP("Momentum") << "  (" << px() << ", " << py() << ", " << pz() << ")\n\t"
                              << "* (" << mom.px() << ", " << mom.py() << ", " << mom.pz() << ")\n\t"
                              << "= " << px() * mom.px() + py() * mom.py() + pz() * mom.pz();
    return px() * mom.px() + py() * mom.py() + pz() * mom.pz();
  }

  double Momentum::fourProduct(const Momentum& mom) const { return (transposed() * metricMatrix * mom)(0); }

  double Momentum::crossProduct(const Momentum& mom) const { return px() * mom.py() - py() * mom.px(); }

  //--- various setters

  Momentum& Momentum::setPx(double px) {
    (*this)(X) = px;
    return computeP();
  }

  Momentum& Momentum::setPy(double py) {
    (*this)(Y) = py;
    return computeP();
  }

  Momentum& Momentum::setPz(double pz) {
    (*this)(Z) = pz;
    return computeP();
  }

  Momentum& Momentum::setEnergy(double e) {
    (*this)(E) = e;
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
    Vector::truncate(tolerance);
    return computeP();
  }

  //--- various getters

  std::array<double, 5> Momentum::pVector() const { return std::array<double, 5>{px(), py(), pz(), energy(), mass()}; }

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

  Momentum Momentum::transverse() const {
    return Momentum::fromPxPyPzE(px(), py(), 0., utils::fastSqrtSqDiff(energy(), pz()));
  }

  double Momentum::eta() const {
    const auto pt_value = pt();
    return (utils::positive(pt_value) ? std::log((p() + fabs(pz())) / pt_value) : 9999.) * utils::sign(pz());
  }

  double Momentum::rapidity() const {
    return utils::positive(energy()) ? std::log((energy() + pz()) / (energy() - pz())) * 0.5
                                     : std::numeric_limits<double>::infinity() * utils::sign(pz());
  }

  double Momentum::deltaEta(const Momentum& oth) const { return fabs(eta() - oth.eta()); }

  double Momentum::deltaPhi(const Momentum& oth) const {
    return normalisePhi(phi() - oth.phi(), {-M_PI, M_PI});  // has to be contained in [-M_PI, M_PI]
  }

  double Momentum::deltaPt(const Momentum& oth) const { return fabs(pt() - oth.pt()); }

  double Momentum::deltaR(const Momentum& oth) const {
    return utils::fastHypot(rapidity() - oth.rapidity(), deltaPhi(oth));
  }

  //--- boosts/rotations

  double Momentum::beta() const {
    const auto mom = p(), ene = energy();
    if (ene == 0.) {
      if (mom == 0.)
        return 0.;
      else {
        CG_WARNING("Momentum:beta") << "beta computed for t=0 momentum.";
        return 1. / ene;
      }
    }
    if (mass2() <= 0.)
      CG_WARNING("Momentum:beta") << "beta computed for an invalid, non-time-like momentum.";
    return mom / ene;
  }

  double Momentum::gamma2() const {
    const auto mom2 = p2(), ene2 = energy2();
    if (ene2 == 0.) {
      if (mom2 == 0.)
        return 1.;
      CG_WARNING("Momentum:gamma") << "gamma computed for t=0 momentum.";
    }
    if (ene2 < mom2) {
      CG_WARNING("Momentum:gamma") << "gamma computed for an invalid space-like momentum.";
      return 0.;
    } else if (ene2 == mom2)
      CG_WARNING("Momentum:gamma") << "gamma computed for a light-like momentum.";
    return ene2 / (ene2 - mom2);
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
    const double pf4 = dot(mom) / mass;
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
    const Matrix rot{
        {+cos_phi, sign * sin_phi, 0., 0.},  // px
        {-sin_phi, sign * cos_phi, 0., 0.},  // py
        {0., 0., 1., 0.},                    // pz
        {0., 0., 0., 1.}                     // e
    };
    *this = rot * (*this);
    return *this;
  }

  Momentum& Momentum::rotateThetaPhi(double theta, double phi) {
    const double cos_theta = std::cos(theta), sin_theta = std::sin(theta);
    const double cos_phi = std::cos(phi), sin_phi = std::sin(phi);
    const Matrix rot{
        {-sin_phi, -cos_theta * cos_phi, sin_theta * cos_phi, 0.},  // px
        {+cos_phi, -cos_theta * sin_phi, sin_theta * sin_phi, 0.},  // py
        {0., sin_theta, cos_theta, 0.},                             // pz
        {0., 0., 0., 1.}                                            // e
    };
    //FIXME check this! cos(phi)->-sin(phi) & sin(phi)->cos(phi) --> phi->phi+pi/2 ?
    *this = rot * (*this);
    return *this;
  }

  //--- printout

  std::ostream& operator<<(std::ostream& os, const Momentum& mom) {
    return os << utils::format(
               "(%8.3f,%8.3f,%8.3f;%8.3f|%8.3f)", mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass());
  }
}  // namespace cepgen
