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
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  /// Express an angle in between two extrema
  static double normalisePhi(double phi, const Limits& range) {
    static const double M_2PI = 2. * M_PI;
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

  static double fastHypot(double x, double y) {
    x *= x, y *= y;
    return std::sqrt(x + y);
  }

  static double fastHypot(double x, double y, double z) {
    x *= x, y *= y, z *= z;
    return std::sqrt(x + y + z);
  }

  Momentum::Momentum(double x, double y, double z, double t) : std::array<double, 4>{{x, y, z, t == -1. ? 0. : t}} {
    computeP();
  }

  Momentum::Momentum(double* p) {
    std::copy(p, p + 4, begin());
    computeP();
  }

  //--- static constructors

  Momentum Momentum::fromPtEtaPhiE(double pt, double eta, double phi, double e) {
    return Momentum(pt * cos(phi), pt * sin(phi), pt * sinh(eta), e);
  }

  Momentum Momentum::fromPtEtaPhiM(double pt, double eta, double phi, double m) {
    return Momentum::fromPtEtaPhiE(pt, eta, phi, fastHypot(pt * cosh(eta), m));
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
    const double et = fastHypot(px, py, m);
    return Momentum(px, py, et * sinh(rap), et * cosh(rap));
  }

  Momentum Momentum::fromPtYPhiM(double pt, double rap, double phi, double m) {
    const double et = fastHypot(pt, m);
    return Momentum(pt * cos(phi), pt * sin(phi), et * sinh(rap), et * cosh(rap));
  }

  //--- arithmetic operators

  Momentum Momentum::operator+(const Momentum& mom) const {
    return Momentum(px() + mom.px(), py() + mom.py(), pz() + mom.pz(), energy() + mom.energy());
  }

  Momentum& Momentum::operator+=(const Momentum& mom) {
    *this = *this + mom;
    return computeP();
  }

  Momentum Momentum::operator-() const { return Momentum(-px(), -py(), -pz(), energy()); }

  Momentum Momentum::operator-(const Momentum& mom) const {
    return Momentum(px() - mom.px(), py() - mom.py(), pz() - mom.pz(), energy() - mom.energy());
  }

  Momentum& Momentum::operator-=(const Momentum& mom) {
    *this = *this - mom;
    return computeP();
  }

  double Momentum::operator*(const Momentum& mom) const { return threeProduct(mom); }

  Momentum Momentum::operator%(const Momentum& mom) const {
    return Momentum(
        py() * mom.pz() - pz() * mom.py(), pz() * mom.px() - px() * mom.pz(), px() * mom.py() - py() * mom.px());
  }

  Momentum Momentum::operator*(double c) const { return Momentum(c * px(), c * py(), c * pz(), c * energy()); }

  Momentum& Momentum::operator*=(double c) {
    *this = *this * c;
    return computeP();
  }

  Momentum operator*(double c, const Momentum& mom) {
    return Momentum(c * mom.px(), c * mom.py(), c * mom.pz(), c * mom.energy());
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

  double Momentum::crossProduct(const Momentum& mom) const { return px() * mom.py() - py() * mom.px(); }

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

  Momentum& Momentum::setMass(double m) { return setEnergy(fastHypot(p_, m)).computeP(); }

  Momentum& Momentum::setMass2(double m2) { return setEnergy(std::sqrt(p_ * p_ + m2)).computeP(); }

  Momentum& Momentum::setP(double px, double py, double pz, double e) { return setEnergy(e).setP(px, py, pz); }

  Momentum& Momentum::setP(double px, double py, double pz) { return setPx(px).setPy(py).setPz(pz); }

  //--- kinematics constrainers

  Momentum& Momentum::computeEnergyFromMass(double on_shell_mass) { return setMass(on_shell_mass); }

  Momentum& Momentum::computePzFromMass(double on_shell_mass) {
    return setPz(normaliseSqrt(pz() * pz() + mass2() - on_shell_mass * on_shell_mass));
  }

  Momentum& Momentum::computeP() {
    p_ = fastHypot(px(), py(), pz());
    return *this;
  }

  Momentum& Momentum::truncate(double tolerance) {
    std::replace_if(
        begin(), end(), [&tolerance](const auto& p) { return p <= tolerance; }, 0.);
    return computeP();
  }

  //--- various getters

  std::array<double, 5> Momentum::pVector() const {
    std::array<double, 5> out;
    std::copy(begin(), end(), out.begin());
    out[4] = mass();
    return out;
  }

  double Momentum::energyT2() const {
    const auto ptsq = pt2();
    return ptsq > 0. ? energy2() * ptsq / (ptsq + pz() * pz()) : 0.;
  }

  double Momentum::energyT() const { return normaliseSqrt(energyT2()); }

  double Momentum::mass2() const { return energy2() - p2(); }

  double Momentum::mass() const { return normaliseSqrt(mass2()); }

  double Momentum::massT2() const { return energy2() - pz() * pz(); }

  double Momentum::massT() const { return normaliseSqrt(massT2()); }

  double Momentum::theta() const { return atan2(pt(), pz()); }

  double Momentum::phi() const { return normalisePhi(atan2(py(), px()), {0., 2. * M_PI}); }

  double Momentum::pt() const { return fastHypot(px(), py()); }

  double Momentum::pt2() const { return px() * px() + py() * py(); }

  Momentum Momentum::transverse() const { return Momentum::fromPxPyPzE(px(), py(), 0., energy()); }

  double Momentum::eta() const {
    const int sign = pz() / fabs(pz());
    const auto ptval = pt();
    return (ptval != 0. ? std::log((p() + fabs(pz())) / ptval) : 9999.) * sign;
  }

  double Momentum::rapidity() const {
    const int sign = pz() / fabs(pz());
    return energy() >= 0. ? std::log((energy() + pz()) / (energy() - pz())) * 0.5 : 999. * sign;
  }

  double Momentum::deltaEta(const Momentum& oth) const { return fabs(eta() - oth.eta()); }

  double Momentum::deltaPhi(const Momentum& oth) const {
    return normalisePhi(phi() - oth.phi(), {-M_PI, M_PI});  // has to be contained in [-M_PI, M_PI]
  }

  double Momentum::deltaPt(const Momentum& oth) const { return fabs(pt() - oth.pt()); }

  double Momentum::deltaR(const Momentum& oth) const { return fastHypot(rapidity() - oth.rapidity(), deltaPhi(oth)); }

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
      CG_WARNING("Momentum:beta") << "beta computed for an invalid, non-timelike momentum.";
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
      CG_WARNING("Momentum:gamma") << "gamma computed for an invalid spacelike momentum.";
      return 0.;
    } else if (ene2 == mom2)
      CG_WARNING("Momentum:gamma") << "gamma computed for a lightlike momentum.";
    return ene2 / (ene2 - mom2);
  }

  double Momentum::gamma() const { return std::sqrt(gamma2()); }

  Momentum& Momentum::betaGammaBoost(double gamma, double betagamma) {
    if (gamma == 1. && betagamma == 0.)
      return *this;  // trivial case

    const double apz = pz(), ae = energy();

    return setEnergy(gamma * ae + betagamma * apz).setPz(gamma * apz + betagamma * ae);
  }

  Momentum& Momentum::lorentzBoost(const Momentum& mom) {
    //--- do not boost on a system at rest
    if (mom.p() == 0.)
      return *this;

    const double mass = mom.mass();
    const double pf4 = ((*this)[X] * mom[X] + (*this)[Y] * mom[Y] + (*this)[Z] * mom[Z] + (*this)[E] * mom[E]) / mass;
    const double fn = (pf4 + (*this)[E]) / (mom[E] + mass);
    (*this) += fn * mom;
    return setEnergy(pf4);
  }

  Momentum& Momentum::rotatePhi(double phi, double sign) {
    const double sphi = sin(phi), cphi = cos(phi);
    const double px = (*this)[X] * cphi + sign * (*this)[Y] * sphi, py = -(*this)[X] * sphi + sign * (*this)[Y] * cphi;
    return setPx(px).setPy(py);
  }

  Momentum& Momentum::rotateThetaPhi(double theta, double phi) {
    const double ctheta = cos(theta), stheta = sin(theta);
    const double cphi = cos(phi), sphi = sin(phi);
    double rotmtx[3][3], mom[3];  //FIXME check this! cos(phi)->-sin(phi) & sin(phi)->cos(phi) --> phi->phi+pi/2 ?
    rotmtx[X][X] = -sphi;
    rotmtx[X][Y] = -ctheta * cphi;
    rotmtx[X][Z] = stheta * cphi;
    rotmtx[Y][X] = cphi;
    rotmtx[Y][Y] = -ctheta * sphi;
    rotmtx[Y][Z] = stheta * sphi;
    rotmtx[Z][X] = 0.;
    rotmtx[Z][Y] = stheta;
    rotmtx[Z][Z] = ctheta;

    for (size_t i = X; i <= Z; ++i) {
      mom[i] = 0.;
      for (size_t j = X; j <= Z; ++j)
        mom[i] += rotmtx[i][j] * (*this)[j];
    }
    return setP(mom[X], mom[Y], mom[Z]);
  }

  //--- printout

  std::ostream& operator<<(std::ostream& os, const Momentum& mom) {
    return os << utils::format(
               "(%8.3f,%8.3f,%8.3f;%8.3f|%8.3f)", mom.px(), mom.py(), mom.pz(), mom.energy(), mom.mass());
  }
}  // namespace cepgen
