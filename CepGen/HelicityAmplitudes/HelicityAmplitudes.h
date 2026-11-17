/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2026  Laurent Forthomme
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

#ifndef CepGen_HelicityAmplitudes_HelicityAmplitude_h
#define CepGen_HelicityAmplitudes_HelicityAmplitude_h

#include <complex>

#include "CepGen/Modules/ModuleFactory.h"
#include "CepGen/Modules/NamedModule.h"

/// Add an amplitudes computation algorithm definition to the factory
#define REGISTER_HELICITY_AMPLITUDES(name, obj)                                           \
  namespace cepgen::amplitudes {                                                          \
    struct BUILDER_NAME(obj) {                                                            \
      BUILDER_NAME(obj)() { HelicityAmplitudesFactory::get().registerModule<obj>(name); } \
    };                                                                                    \
    static const BUILDER_NAME(obj) gNachAmp##obj;                                         \
  }                                                                                       \
  static_assert(true, "")

namespace cepgen {
  class ParametersList;
}  // namespace cepgen

namespace cepgen::amplitudes {
  /// Helper container to handle all kinematics variables computation once
  class Kinematics {
  public:
    explicit Kinematics(double mv2, double s_hat, double t_hat, double u_hat);
    static Kinematics fromSCosTheta(double shat, double cos_theta, double mv2);

    bool operator!=(const Kinematics&) const;
    friend std::ostream& operator<<(std::ostream&, const Kinematics&);

    // base variables
    const double shat{0.}, that{0.}, uhat{0.};

  private:
    void setCosTheta(double cos_theta);
    const double mv2_{0.};  /// Squared mass of each of the two produced bosons (in GeV^2)

  public:
    // all derived variables
    const double shat2{0.};
    const double beta2{0.}, beta{0.};
    const double inv_gamma2{0.}, gamma2{0.}, gamma{0.}, inv_gamma{0.};
    double cos_theta{0.}, cos_theta2{0.}, sin_theta2{0.}, sin_theta{0.};
    double invA{0.};
  };

  /// Simple container for helicity components
  struct HelicityStates {
    short lam1;  ///< first incoming photon
    short lam2;  ///< second incoming photon
    short lam3;  ///< first outgoing boson
    short lam4;  ///< second outgoing boson
  };

  /// Amplitudes computational tool, as developed by Nachtmann et al.
  /// \cite Nachtmann:2005en
  class HelicityAmplitudes : public NamedModule<HelicityAmplitudes> {
  public:
    static ParametersDescription description();

    /// Compute the amplitude for a given kinematics and a given set of helicity components
    std::complex<double> operator()(const Kinematics&, short lam1, short lam2, short lam3, short lam4) const;

    virtual bool hasAAWW() const { return false; }
    virtual bool hasAAZZ() const { return false; }

  protected:
    explicit HelicityAmplitudes(const ParametersList&);

    virtual std::complex<double> amplitude(const Kinematics&, const HelicityStates&) const = 0;
    double s1() const { return s1_; }
    double c1() const { return std::sqrt(1. - s1_ * s1_); }

  private:
    const double s1_;
  };

  /// A helicity-dependent amplitudes calculation algorithms factory
  DEFINE_FACTORY(HelicityAmplitudesFactory, HelicityAmplitudes, "Helicity-dependent amplitudes factory");
}  // namespace cepgen::amplitudes

#endif
