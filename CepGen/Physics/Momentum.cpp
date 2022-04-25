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

namespace cepgen {
  Momentum::Momentum(double x, double y, double z, double t) : std::array<double, 4>{{x, y, z, t == -1. ? 0. : t}} {
    computeP();
  }

  Momentum::Momentum(double* p) {
    std::copy(p, p + 4, begin());
    computeP();
  }

  //--- static constructors

  Momentum Momentum::fromPtEtaPhiE(double pt, double eta, double phi, double e) {
    const double px = pt * cos(phi), py = pt * sin(phi), pz = pt * sinh(eta);
    return Momentum(px, py, pz, e);
  }

  Momentum Momentum::fromPtEtaPhiM(double pt, double eta, double phi, double m) {
    const double px = pt * cos(phi), py = pt * sin(phi), pz = pt * sinh(eta), e = px * px + py * py + pz * pz + m * m;
    return Momentum(px, py, pz, e);
  }

  Momentum Momentum::fromPThetaPhiE(double p, double theta, double phi, double e) {
    const double px = p * sin(theta) * cos(phi), py = p * sin(theta) * sin(phi), pz = p * cos(theta);
    return Momentum(px, py, pz, e);
  }

  Momentum Momentum::fromPxPyPzE(double px, double py, double pz, double e) { return Momentum(px, py, pz, e); }

  Momentum Momentum::fromPxPyPzM(double px, double py, double pz, double m) {
    return Momentum(px, py, pz).setMass(m).computeP();
  }

  Momentum Momentum::fromPxPyYM(double px, double py, double rap, double m) {
    const double pt = std::hypot(px, py), et = std::hypot(pt, m);
    return Momentum(px, py, et * sinh(rap), et * cosh(rap));
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

  Momentum& Momentum::setMass2(double m2) { return setEnergy(sqrt(p2() + m2)).computeP(); }

  Momentum& Momentum::setP(double px, double py, double pz, double e) { return setP(px, py, pz).setEnergy(e); }

  Momentum& Momentum::setP(double px, double py, double pz) { return setPx(px).setPy(py).setPz(pz).computeP(); }

  Momentum& Momentum::computeP() {
    p_ = std::hypot(pt(), pz());
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

  double Momentum::energyT() const {
    const auto et2 = energyT2();
    return et2 > 0. ? std::sqrt(et2) : -std::sqrt(-et2);
  }

  double Momentum::mass() const {
    if (mass2() >= 0.)
      return sqrt(mass2());
    return -sqrt(-mass2());
  }

  double Momentum::massT2() const { return energy2() - pz() * pz(); }

  double Momentum::massT() const {
    const auto mt2 = massT2();
    return mt2 > 0. ? std::sqrt(mt2) : -std::sqrt(-mt2);
  }

  double Momentum::theta() const { return atan2(pt(), pz()); }

  double Momentum::phi() const { return atan2(py(), px()); }

  double Momentum::pt() const { return std::hypot(px(), py()); }

  double Momentum::pt2() const { return px() * px() + py() * py(); }

  double Momentum::eta() const {
    const int sign = (pz() / fabs(pz()));
    return pt() != 0. ? log((p() + fabs(pz())) / pt()) * sign : 9999. * sign;
  }

  double Momentum::rapidity() const {
    const int sign = (pz() / fabs(pz()));
    return energy() >= 0. ? log((energy() + pz()) / (energy() - pz())) * 0.5 : 999. * sign;
  }

  double Momentum::deltaEta(const Momentum& oth) const { return fabs(eta() - oth.eta()); }

  double Momentum::deltaPhi(const Momentum& oth) const {
    static const double M_2PI = 2. * M_PI;
    double dphi = phi() - oth.phi();
    // has to be contained in [-M_PI, M_PI]
    while (dphi < -M_PI)
      dphi += M_2PI;
    while (dphi > +M_PI)
      dphi -= M_2PI;
    return dphi;
  }

  double Momentum::deltaPt(const Momentum& oth) const { return fabs(pt() - oth.pt()); }

  double Momentum::deltaR(const Momentum& oth) const { return std::hypot(rapidity() - oth.rapidity(), deltaPhi(oth)); }

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

  double Momentum::gamma() const {
    const auto mom2 = p2(), ene2 = energy2();
    if (ene2 == 0.) {
      if (mom2 == 0.) {
        return 1.;
      } else
        CG_WARNING("Momentum:gamma") << "gamma computed for t=0 momentum.";
    }
    if (ene2 < mom2) {
      CG_WARNING("Momentum:gamma") << "gamma computed for an invalid spacelike momentum.";
      return 0.;
    } else if (ene2 == mom2)
      CG_WARNING("Momentum:gamma") << "gamma computed for a lightlike momentum.";
    return 1. / std::sqrt(1. - mom2 / ene2);
  }

  Momentum& Momentum::betaGammaBoost(double gamma, double betagamma) {
    if (gamma == 1. && betagamma == 0.)
      return *this;  // trivial case

    const double apz = pz(), ae = energy();

    return setPz(gamma * apz + betagamma * ae).setEnergy(gamma * ae + betagamma * apz).computeP();
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
    return os << "(" << mom.energy() << "|" << mom.px() << " " << mom.py() << " " << mom.pz() << ")";
  }

  double CMEnergy(const Momentum& m1, const Momentum& m2) {
    if (m1.mass() * m2.mass() < 0. || m1.energy() * m2.energy() < 0.)
      return 0.;
    return sqrt(m1.mass2() + m2.mass2() + 2. * m1.energy() * m2.energy() - 2. * (m1 * m2));
  }
}  // namespace cepgen
