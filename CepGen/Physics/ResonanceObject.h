/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGen_Physics_ResonanceObject_h
#define CepGen_Physics_ResonanceObject_h

#include <vector>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  /// General definition for a resonance
  class ResonanceObject : public SteeredObject<ResonanceObject> {
  public:
    explicit ResonanceObject(const ParametersList&);

    static ParametersDescription description();

    /// kinematics needed for threshold relativistic B-W
    struct KinematicsBlock {
      explicit KinematicsBlock(double w2, double q2, double mp2, double mpi2, double meta2);
      static double mom(double energy, double mass2) { return std::sqrt(std::max(0., energy * energy - mass2)); }
      const double w2, w;
      const double q2;
      // equivalent photon energy-momentum
      const double k, kcm;
      // pion momentum
      const double ppicm;
      // two-pion momentum
      const double ppi2cm;
      // eta meson momentum
      const double petacm;
    };

  protected:
    double kr() const { return 0.5 * (mass_ * mass_ - mp2_) / mp_; }
    double pcmr(double m2) const { return KinematicsBlock::mom(ecmr(m2), m2); }
    double ecmr(double m2) const;
    double kcmr() const { return ecmr(0.); }
    /// partial widths for all decays
    double partialWidth(const KinematicsBlock&) const;
    /// virtual photon width
    double photonWidth(const KinematicsBlock&) const;

    /// Branching ratios container for resonance decay into single, double pion or eta states
    const struct BranchingRatios : SteeredObject<BranchingRatios> {
      explicit BranchingRatios(const ParametersList&);

      static ParametersDescription description();

      bool valid() const { return singlepi + doublepi + eta == 1.; }
      /// single pion branching ratio
      double singlepi;
      /// double pion branching ratio
      double doublepi;
      /// eta meson branching ratio
      double eta;
    } br_;
    /// meson angular momentum
    const int ang_mom_;
    /// damping parameter
    const double x0_;
    /// mass, in GeV/c2
    const double mass_;
    /// full width, in GeV
    const double width_;
    /// proton mass, in GeV/c^2
    const double mp_;
    /// proton squared mass, in GeV^2/c^4
    const double mp2_;
    /// pion squared mass, in GeV^2/c^4
    const double mpi2_;
    /// eta meson squared mass, in GeV^2/c^4
    const double meta2_;
    /// squared damping parameter
    const double x02_;
  };
}  // namespace cepgen

#endif
