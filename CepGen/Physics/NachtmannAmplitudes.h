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

    /// Helper container to handle all kinematics variables computation once
    class Kinematics {
    public:
      Kinematics(double mw2, double shat, double that, double uhat);
      bool operator!=(const Kinematics&) const;

      // base variables
      const double shat, that, uhat;

    private:
      /// W squared mass, in GeV^2
      const double mw2_;

    public:
      // all derived variables
      const double beta2, beta;
      const double inv_gamma2, gamma2, gamma, inv_gamma;
      const double cos_theta, cos_theta2, sin_theta2, sin_theta;
      const double invA;
    };

    /// Compute the amplitude for a given kinematics and a given set of helicity components
    std::complex<double> operator()(const Kinematics&, short lam1, short lam2, short lam3, short lam4) const;

  private:
    const Mode mode_;
    /// Collection of parameters for the EFT extension
    const struct EFTParameters {
      explicit EFTParameters(const ParametersList& params);
      const double s1, mH;
      double c1() const { return sqrt(1. - s1 * s1); }
    } eft_ext_;

    /// Simple container for helicity components
    struct Helicities {
      short lam1;  ///< first incoming photon
      short lam2;  ///< second incoming photon
      short lam3;  ///< first outgoing W
      short lam4;  ///< second outgoing W
    };
    const double G_EM_SQ;
    const double G_EM;

    /// Compute the amplitude for the Standard model
    std::complex<double> amplitudeSM(const Kinematics&, const Helicities&) const;
    std::complex<double> amplitudeW(const Kinematics&, const Helicities&) const;
    std::complex<double> amplitudeWbar(const Kinematics&, const Helicities&) const;
    std::complex<double> amplitudephiW(const Kinematics&, const Helicities&) const;
    std::complex<double> amplitudeWB(const Kinematics&, const Helicities&) const;
    std::complex<double> amplitudeWbarB(const Kinematics&, const Helicities&) const;
  };
}  // namespace cepgen
