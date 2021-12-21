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

#ifndef CepGen_Physics_Momentum_h
#define CepGen_Physics_Momentum_h

#include <array>
#include <iosfwd>

namespace cepgen {
  /**
   * Container for a particle's 4-momentum, along with useful methods to ease the development of any matrix element level generator
   * \brief 4-momentum for a particle
   * \date Dec 2015
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   */
  class Momentum : private std::array<double, 4> {
  public:
    /// Build a 4-momentum at rest with an invalid energy (no mass information known)
    Momentum();
    /// Build a 4-momentum using its 3-momentum coordinates and its energy
    Momentum(double x, double y, double z, double t = -1.);
    /// Build a 4-momentum using its 3-momentum coordinates and its energy
    Momentum(double* p);

    //--- static definitions

    /// Build a 3-momentum from its three pseudo-cylindrical coordinates
    static Momentum fromPtEtaPhiE(double pt, double eta, double phi, double e = -1.);
    /// Build a 3-momentum from its three pseudo-cylindrical coordinates
    static Momentum fromPtEtaPhiM(double pt, double eta, double phi, double m);
    /// Build a 4-momentum from its scalar momentum, and its polar and azimuthal angles
    static Momentum fromPThetaPhiE(double p, double theta, double phi, double e = -1.);
    /// Build a 4-momentum from its four momentum and energy coordinates
    static Momentum fromPxPyPzE(double px, double py, double pz, double e);
    /// Build a 4-momentum from its three momentum coordinates and mass
    static Momentum fromPxPyPzM(double px, double py, double pz, double m);
    /// Build a 4-momentum from its transverse momentum, rapidity and mass
    static Momentum fromPxPyYM(double px, double py, double rap, double m);

    //--- vector and scalar operators

    /// Scalar product of the 3-momentum with another 3-momentum
    double threeProduct(const Momentum&) const;
    /// Scalar product of the 4-momentum with another 4-momentum
    double fourProduct(const Momentum&) const;
    /// Vector product of the 3-momentum with another 3-momentum
    double crossProduct(const Momentum&) const;
    /// Compute the 4-vector sum of two 4-momenta
    Momentum operator+(const Momentum&) const;
    /// Add a 4-momentum through a 4-vector sum
    Momentum& operator+=(const Momentum&);
    /// Unary inverse operator
    Momentum operator-() const;
    /// Compute the inverse per-coordinate 4-vector
    Momentum operator-(const Momentum&) const;
    /// Subtract a 4-momentum through a 4-vector sum
    Momentum& operator-=(const Momentum&);
    /// Scalar product of two 3-momenta
    double operator*(const Momentum&) const;
    /// Vector product of two 3-momenta
    Momentum operator%(const Momentum&) const;
    /// Scalar product of the 3-momentum with another 3-momentum
    double operator*=(const Momentum&);
    /// Multiply all components of a 4-momentum by a scalar
    Momentum operator*(double c) const;
    /// Multiply all 4-momentum coordinates by a scalar
    Momentum& operator*=(double c);
    /// Left-multiply all 4-momentum coordinates by a scalar
    friend Momentum operator*(double, const Momentum&);
    /// Human-readable format for a particle's momentum
    friend std::ostream& operator<<(std::ostream&, const Momentum&);

    /// Forward \f$\beta-\gamma\f$ boost
    Momentum& betaGammaBoost(double gamma, double betagamma);
    /// Forward Lorentz boost
    Momentum& lorentzBoost(const Momentum& p);

    //--- setters and getters

    /// Set all the components of the 4-momentum (in GeV)
    Momentum& setP(double px, double py, double pz, double e);
    /// Set all the components of the 3-momentum (in GeV)
    Momentum& setP(double px, double py, double pz);
    /// Set the momentum along the \f$x\f$-axis (in GeV)
    inline Momentum& setPx(double px) {
      (*this)[X] = px;
      return *this;
    }
    /// Momentum along the \f$x\f$-axis (in GeV)
    inline double px() const { return (*this)[X]; }
    /// Set the momentum along the \f$y\f$-axis (in GeV)
    inline Momentum& setPy(double py) {
      (*this)[Y] = py;
      return *this;
    }
    /// Momentum along the \f$y\f$-axis (in GeV)
    inline double py() const { return (*this)[Y]; }
    /// Set the longitudinal momentum (in GeV)
    inline Momentum& setPz(double pz) {
      (*this)[Z] = pz;
      return *this;
    }
    /// Longitudinal momentum (in GeV)
    inline double pz() const { return (*this)[Z]; }
    /// Transverse momentum (in GeV)
    double pt() const;
    /// Squared transverse momentum (in GeV\f$^2\f$)
    double pt2() const;
    /// 5-vector of double precision floats (in GeV)
    std::array<double, 5> pVector() const;
    /// 3-momentum norm (in GeV)
    inline double p() const { return p_; }
    /// Squared 3-momentum norm (in GeV\f$^2\f$)
    inline double p2() const { return p_ * p_; }
    /// Set the energy (in GeV)
    inline Momentum& setEnergy(double e) {
      (*this)[E] = e;
      return *this;
    }
    /// Energy (in GeV)
    inline double energy() const { return (*this)[E]; }
    /// Squared energy (in GeV\f$^2\f$)
    inline double energy2() const { return (*this)[E] * (*this)[E]; }
    /// Compute the energy from the mass
    Momentum& setMass2(double m2);
    /// Squared mass (in GeV\f$^2\f$) as computed from its energy and momentum
    inline double mass2() const { return energy2() - p2(); }
    /// Compute the energy from the mass
    inline Momentum& setMass(double m) { return setMass2(m * m); }
    /// Mass (in GeV) as computed from its energy and momentum
    /// \note Returns \f$-\sqrt{|E^2-\mathbf{p}^2|}<0\f$ if \f$\mathbf{p}^2>E^2\f$
    double mass() const;
    /// Polar angle (angle with respect to the longitudinal direction)
    double theta() const;
    /// Azimuthal angle (angle in the transverse plane)
    double phi() const;
    /// Pseudo-rapidity
    double eta() const;
    /// Rapidity
    double rapidity() const;

    /// Beta scalar value
    double beta() const;
    /// Gamma scalar value
    double gamma() const;

    /// Apply a threshold to all values with a given tolerance
    Momentum& truncate(double tolerance = 1.e-10);
    /// Rotate the transverse components by an angle phi (and reflect the y coordinate)
    Momentum& rotatePhi(double phi, double sign);
    /// Rotate the particle's momentum by a polar/azimuthal angle
    Momentum& rotateThetaPhi(double theta_, double phi_);
    /// Apply a \f$ x\rightarrow -x\f$ transformation
    inline Momentum& mirrorX() {
      (*this)[X] *= -1.;
      return *this;
    }
    /// Apply a \f$ y\rightarrow -y\f$ transformation
    inline Momentum& mirrorY() {
      (*this)[Y] *= -1.;
      return *this;
    }
    /// Apply a \f$ z\rightarrow -z\f$ transformation
    inline Momentum& mirrorZ() {
      (*this)[Z] *= -1.;
      return *this;
    }

  private:
    enum coord_t { X = 0, Y = 1, Z = 2, E = 3 };
    /// Compute the 3-momentum's norm
    Momentum& computeP();
    /// 3-momentum's norm (in GeV/c)
    double p_{0.};
  };
  /// Compute the centre of mass energy of two particles momenta
  double CMEnergy(const Momentum& m1, const Momentum& m2);
}  // namespace cepgen

#endif
