/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#ifndef CepGen_Physics_NachtmannAmplitudes_h
#define CepGen_Physics_NachtmannAmplitudes_h

#include <complex>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  class ParametersList;
  /// Amplitudes computational tool, as developed by Nachtmann et al.
  /// \cite Nachtmann:2005en
  class NachtmannAmplitudes final : public SteeredObject<NachtmannAmplitudes> {
  public:
    explicit NachtmannAmplitudes(const ParametersList&);

    /// Model giving an amplitude for the two-photon WW production
    enum class Mode { SM, W, Wbar, phiW, phiWbar, phiB, phiBbar, WB, WbarB };
    friend std::ostream& operator<<(std::ostream&, const Mode&);
    const Mode& mode() const { return mode_; }

    /// Helper container to handle all kinematics variables computation once
    class Kinematics {
    public:
      explicit Kinematics(double mw2, double shat, double that, double uhat);
      bool operator!=(const Kinematics&) const;
      friend std::ostream& operator<<(std::ostream&, const Kinematics&);
      static Kinematics fromScosTheta(double shat, double cos_theta, double mw2);

      // base variables
      const double shat{0.}, that{0.}, uhat{0.};

    private:
      void setCosTheta(double cos_theta);
      /// W squared mass, in GeV^2
      const double mw2_{0.};

    public:
      // all derived variables
      const double shat2{0.};
      const double beta2{0.}, beta{0.};
      const double inv_gamma2{0.}, gamma2{0.}, gamma{0.}, inv_gamma{0.};
      double cos_theta{0.}, cos_theta2{0.}, sin_theta2{0.}, sin_theta{0.};
      double invA{0.};
    };

    /// Compute the amplitude for a given kinematics and a given set of helicity components
    std::complex<double> operator()(const Kinematics&, short lam1, short lam2, short lam3, short lam4) const;

    static ParametersDescription description();

  private:
    const Mode mode_;
    /// Collection of parameters for the EFT extension
    const struct EFTParameters final : SteeredObject<EFTParameters> {
      explicit EFTParameters(const ParametersList&);
      const double s1, mH;
      double c1() const { return sqrt(1. - s1 * s1); }
      static ParametersDescription description();
    } eft_ext_;
    /// Simple container for helicity components
    struct HelicityStates {
      short lam1;  ///< first incoming photon
      short lam2;  ///< second incoming photon
      short lam3;  ///< first outgoing W
      short lam4;  ///< second outgoing W
    };
    const double G_EM_SQ;
    const double G_EM;

    /// Compute the amplitude for the Standard model
    std::complex<double> amplitudeSM(const Kinematics&, const HelicityStates&) const;
    std::complex<double> amplitudeW(const Kinematics&, const HelicityStates&) const;
    std::complex<double> amplitudeWbar(const Kinematics&, const HelicityStates&) const;
    std::complex<double> amplitudephiW(const Kinematics&, const HelicityStates&) const;
    std::complex<double> amplitudeWB(const Kinematics&, const HelicityStates&) const;
    std::complex<double> amplitudeWbarB(const Kinematics&, const HelicityStates&) const;
  };
}  // namespace cepgen

#endif
