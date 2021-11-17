/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *                2017-2019  Wolfgang Schaefer
 *                2019       Marta Luszczak
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

#include <complex>

namespace cepgen {
  class ParametersList;
  class NachtmannAmplitudes {
  public:
    NachtmannAmplitudes(const ParametersList&);

    /// Model giving an amplitude for the two-photon WW production
    enum class Mode { SM, W, Wbar, phiW, phiWbar, phiB, phiBbar, WB, WbarB };
    friend std::ostream& operator<<(std::ostream&, const Mode&);
    const Mode& mode() const { return mode_; }

    /// Compute the amplitude for a given kinematics and a given set of helicity components
    std::complex<double> operator()(
        double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4) const;

  private:
    const Mode mode_;
    /// Collection of parameters for the EFT extension
    const struct EFTParameters {
      explicit EFTParameters(const ParametersList& params);
      const double s1, mH;
      double c1() const { return sqrt(1. - s1 * s1); }
    } eft_ext_;

    /// Helper container to handle all kinematics variables computation once
    class Kinematics {
    public:
      Kinematics(double mw2, double shat, double that, double uhat);
      bool isEqual(double shat, double that, double uhat) const;

      // base variables
      double shat, that, uhat;

    private:
      double mw2_;

    public:
      // all derived variables
      double beta2, beta;
      double inv_gamma2, gamma2, gamma, inv_gamma;
      double cos_theta, cos_theta2, sin_theta2, sin_theta;
      double invA;
    };

    /// Simple container for helicity components
    struct Helicities {
      short lam1;  ///< first incoming photon
      short lam2;  ///< second incoming photon
      short lam3;  ///< first outgoing W
      short lam4;  ///< second outgoing W
    };

    /// Compute the amplitude for the Standard model
    std::complex<double> amplitudeSM(const Helicities&) const;
    std::complex<double> amplitudeW(const Helicities&) const;
    std::complex<double> amplitudeWbar(const Helicities&) const;
    std::complex<double> amplitudephiW(const Helicities&) const;
    std::complex<double> amplitudeWB(const Helicities&) const;
    std::complex<double> amplitudeWbarB(const Helicities&) const;

    /// W squared mass, in GeV^2
    const double mw2_;
    mutable Kinematics kin_{mw2_, 0., 0., 0.};
  };
}  // namespace cepgen
