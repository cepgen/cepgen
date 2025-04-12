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

#ifndef CepGen_Physics_Momentum_h
#define CepGen_Physics_Momentum_h

#include <array>
#include <iosfwd>

namespace cepgen {
  class Vector;
  /// Container for a particle's 4-momentum, along with useful methods to ease the development of any matrix element level generator
  /// \date Dec 2015
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  class Momentum : std::array<double, 4> {
  public:
    explicit Momentum(double x = 0.,
                      double y = 0.,
                      double z = 0.,
                      double t = -1.);  ///< Build a 4-momentum using its 3-momentum coordinates and its energy
    explicit Momentum(double* p);       ///< Build a 4-momentum using its 3-momentum coordinates and its energy
    explicit Momentum(const Vector&);   ///< Build a 4-momentum using its 3-momentum coordinates and its energy

    bool operator==(const Momentum&) const;                                         ///< Equality operator
    inline bool operator!=(const Momentum& oth) const { return !operator==(oth); }  ///< Inequality operator

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
    /// Build a 4-momentum from its transverse momentum, azimuthal angle, rapidity and mass
    static Momentum fromPtYPhiM(double pt, double rap, double phi, double m);
    /// Build a 4-momentum from its transverse momentum, azimuthal angle, rapidity and energy
    static Momentum fromPtYPhiE(double pt, double rap, double phi, double e);

    //--- vector and scalar operators

    double threeProduct(const Momentum&) const;  ///< Scalar product of the 3-momentum with another 3-momentum
    double fourProduct(const Momentum&) const;   ///< Scalar product of the 4-momentum with another 4-momentum
    double crossProduct(const Momentum&) const;  ///< Vector product of the 3-momentum with another 3-momentum

    Momentum operator+(const Momentum&) const;           ///< Compute the 4-vector sum of two 4-momenta
    Momentum& operator+=(const Momentum&);               ///< Add a 4-momentum through a 4-vector sum
    Momentum operator-() const;                          ///< Unary inverse operator
    Momentum operator-(const Momentum&) const;           ///< Compute the inverse per-coordinate 4-vector
    Momentum& operator-=(const Momentum&);               ///< Subtract a 4-momentum through a 4-vector sum
    double operator*(const Momentum&) const;             ///< Scalar product of two 3-momenta
    Momentum operator%(const Momentum&) const;           ///< Vector product of two 3-momenta
    double operator*=(const Momentum&);                  ///< Scalar product of the 3-momentum with another 3-momentum
    Momentum operator*(double) const;                    ///< Multiply all components of a 4-momentum by a scalar
    Momentum& operator*=(double);                        ///< Multiply all 4-momentum coordinates by a scalar
    friend Momentum operator*(double, const Momentum&);  ///< Left-multiply all 4-momentum coordinates by a scalar
    Momentum operator/(double) const;                    ///< Divide all components of a 4-momentum by a scalar
    Momentum& operator/=(double);                        ///< Divide all 4-momentum coordinates by a scalar

    friend std::ostream& operator<<(std::ostream&, const Momentum&);  ///< Human-readable printout of a momentum

    operator Vector() const;  ///< Cast the 4-momentum object into a 4-dimensional vector

    //--- setters and getters

    Momentum& setP(double px, double py, double pz, double e);  ///< Set all the components of the 4-momentum (in GeV)
    Momentum& setP(double px, double py, double pz);            ///< Set all the components of the 3-momentum (in GeV)
    Momentum& setPx(double px);                                 ///< Set the momentum along the \f$x\f$-axis (in GeV)
    Momentum& setPy(double py);                                 ///< Set the momentum along the \f$y\f$-axis (in GeV)
    Momentum& setPz(double pz);                                 ///< Set the longitudinal momentum (in GeV)
    inline double px() const { return (*this)[X]; }             ///< Momentum along the \f$x\f$-axis (in GeV)
    inline double py() const { return (*this)[Y]; }             ///< Momentum along the \f$y\f$-axis (in GeV)
    inline double pz() const { return (*this)[Z]; }             ///< Longitudinal momentum (in GeV)
    double pt() const;                                          ///< Transverse momentum (in GeV)
    double pt2() const;                                         ///< Squared transverse momentum (in GeV\f$^2\f$)

    std::array<double, 5> pVector() const;                             ///< 5-vector of double precision floats (in GeV)
    inline double p() const { return p_; }                             ///< 3-momentum norm (in GeV)
    inline double p2() const { return p_ * p_; }                       ///< Squared 3-momentum norm (in GeV\f$^2\f$)
    Momentum& setEnergy(double);                                       ///< Set the energy (in GeV)
    inline double energy() const { return (*this)[E]; }                ///< Energy (in GeV)
    inline double energy2() const { return (*this)[E] * (*this)[E]; }  ///< Squared energy (in GeV\f$^2\f$)
    double energyT() const;                                            ///< Transverse energy component (in GeV)
    double energyT2() const;     ///< Squared transverse energy component (in GeV\f$^2\f$)
    Momentum& setMass2(double);  ///< Compute the energy from the mass
    double mass2() const;        ///< Squared mass (in GeV\f$^2\f$) as computed from its energy and momentum
    Momentum& setMass(double);   ///< Compute the energy from the mass
    /// Mass (in GeV) as computed from its energy and momentum
    /// \note Returns \f$-\sqrt{|E^2-\mathbf{p}^2|}<0\f$ if \f$\mathbf{p}^2>E^2\f$
    double mass() const;
    double massT2() const;    ///< Squared transverse mass (in GeV\f$^2\f$)
    double massT() const;     ///< Transverse mass (in GeV)
    double theta() const;     ///< Polar angle (angle with respect to the longitudinal direction)
    double phi() const;       ///< Azimuthal angle (angle in the transverse plane)
    double eta() const;       ///< Pseudo-rapidity
    double rapidity() const;  ///< Rapidity

    enum coord_t { X = 0, Y = 1, Z = 2, E = 3 };  ///< Coordinates names
    Momentum transverse(coord_t = Z) const;       ///< Transverse coordinates of a momentum
    Momentum longitudinal(coord_t = Z) const;     ///< Longitudinal component of a momentum

    Momentum& betaGammaBoost(double gamma, double beta_gamma);  ///< Forward \f$\beta-\gamma\f$ boost
    Momentum& lorentzBoost(const Momentum& p);                  ///< Forward Lorentz boost

    double deltaEta(const Momentum&) const;  ///< Pseudo-rapidity distance between two momenta
    double deltaPhi(const Momentum&) const;  ///< Azimuthal angle opening between two momenta
    double deltaPt(const Momentum&) const;   ///< Transverse momentum distance between two momenta
    double deltaR(const Momentum&) const;    ///< Angular distance between two momenta

    double beta() const;    ///< Beta scalar value
    double gamma2() const;  ///< Squared gamma scalar value
    double gamma() const;   ///< Gamma scalar value

    /// Compute the mass from 4-momentum
    /// \param[in] on_shell_mass Specify on-shell mass (in GeV) to constrain energy
    Momentum& computeEnergyFromMass(double on_shell_mass);
    /// Compute the longitudinal coordinate from energy-mass-transverse momentum constraints
    /// \param[in] on_shell_mass Specify on-shell mass (in GeV) to constrain longitudinal momentum
    Momentum& computePzFromMass(double on_shell_mass);

    Momentum& truncate(double tolerance = 1.e-10);  ///< Apply a threshold to all values with a given tolerance

    Momentum& rotatePhi(double phi, double sign);        ///< Rotate transverse components (+ reflect the y coordinate)
    Momentum& rotateThetaPhi(double theta, double phi);  ///< Rotate the particle's momentum by a polar/azimuthal angle
    /// Apply an \f$ x\rightarrow -x\f$ transformation
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
    Momentum& computeP();  ///< Compute the 3-momentum's norm
    double p_{0.};         ///< 3-momentum's norm (in GeV/c)
  };
}  // namespace cepgen

#endif
