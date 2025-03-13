/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  /// General definition for a resonance
  class ResonanceObject : public SteeredObject<ResonanceObject> {
  public:
    explicit ResonanceObject(const ParametersList&);

    static ParametersDescription description();

    /// Kinematics needed for threshold relativistic B-W
    struct KinematicsBlock {
      explicit KinematicsBlock(double w2, double q2, double mp2, double mpi2, double meta2);
      static inline double mom(double energy, double mass2) { return std::sqrt(std::max(0., energy * energy - mass2)); }
      const double w2, w;
      const double q2;
      // equivalent photon energy-momentum
      const double k, kcm;
      const double ppicm;   ///< pion momentum
      const double ppi2cm;  ///< two-pion momentum
      const double petacm;  ///< eta meson momentum
    };

  protected:
    inline double kr() const { return 0.5 * (mass_ * mass_ - mp2_) / mp_; }
    inline double pcmr(double m2) const { return KinematicsBlock::mom(ecmr(m2), m2); }
    double ecmr(double m2) const;
    inline double kcmr() const { return ecmr(0.); }
    double partialWidth(const KinematicsBlock&) const;  ///< partial widths for all decays
    double photonWidth(const KinematicsBlock&) const;   ///< virtual photon width

    /// Branching ratios container for resonance decay into single, double pion or eta states
    const struct BranchingRatios final : SteeredObject<BranchingRatios> {
      explicit BranchingRatios(const ParametersList&);

      static ParametersDescription description();

      /// Sanity check to ensure only three decay channels are opened
      inline bool valid() const { return single_pion + double_pion + eta == 1.; }

      double single_pion;  ///< single pion branching ratio
      double double_pion;  ///< double pion branching ratio
      double eta;          ///< eta meson branching ratio
    } br_;
    const int ang_mom_;   ///< meson angular momentum
    const double x0_;     ///< damping parameter
    const double mass_;   ///< mass, in GeV/c2
    const double width_;  ///< full width, in GeV
    const double mp_;     ///< proton mass, in GeV/c^2
    const double mp2_;    ///< proton squared mass, in GeV^2/c^4
    const double mpi2_;   ///< pion squared mass, in GeV^2/c^4
    const double meta2_;  ///< eta meson squared mass, in GeV^2/c^4
    const double x02_;    ///< squared damping parameter
  };
}  // namespace cepgen

#endif
