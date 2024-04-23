/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2024  Laurent Forthomme
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

#ifndef CepGen_Physics_Cuts_h
#define CepGen_Physics_Cuts_h

#include "CepGen/Event/Particle.h"

namespace cepgen {
  class Event;
  class ParametersList;

  /// A namespace for all kinematic cuts
  namespace cuts {
    /// Centrally produced particles phase space cuts
    struct Central final : public SteeredObject<Central> {
      Central();
      explicit Central(const ParametersList&);

      static ParametersDescription description();
      bool contain(const Particles&, const Event* evt = nullptr) const;

      Limits pt_single;        ///< single particle transverse momentum
      Limits eta_single;       ///< single particle pseudo-rapidity
      Limits rapidity_single;  ///< single particle rapidity
      Limits energy_single;    ///< single particle energy
      Limits mass_single;      ///< single particle mass
      Limits pt_sum;           ///< multiparticle system transverse momentum
      Limits eta_sum;          ///< multiparticle system pseudo-rapidity
      Limits energy_sum;       ///< multiparticle system energy
      Limits mass_sum;         ///< multiparticle system invariant mass
      Limits pt_diff;          ///< transverse momentum balance between the central particles
      Limits phi_diff;         ///< azimuthal angles difference between the central particles
      Limits rapidity_diff;    ///< rapidity balance between the central particles
    };

    /// Initial parton-like particles phase space cuts
    struct Initial final : public SteeredObject<Initial> {
      explicit Initial(const ParametersList&);

      static ParametersDescription description();
      bool contain(const Particles&, const Event* evt = nullptr) const;

      void setParameters(const ParametersList&) override;

      Limits q2_1;  ///< positive-z parton virtuality
      Limits q2_2;  ///< negative-z parton virtuality
      Limits qt;    ///< parton transverse virtuality
      Limits phi;   ///< parton azimuthal angle
    };

    /// Outgoing beam remnant-like particles phase space cuts
    struct Remnants final : public SteeredObject<Remnants> {
      explicit Remnants(const ParametersList&);

      static ParametersDescription description();
      bool contain(const Particles&, const Event* evt = nullptr) const;

      void setParameters(const ParametersList&) override;

      Limits mx;  ///< diffractive mass
      Limits yj;  ///< diffractive jet rapidity
      Limits xi;  ///< longitudinal momentum fraction

      static constexpr double MX_MIN = 1.07;  ///< Minimal diffractive mass for dissociative proton treatment
    };
  }  // namespace cuts
  /// Collection of cuts to be applied on all particle with a given PDG id
  typedef std::unordered_map<pdgid_t, cuts::Central> PerIdCuts;
  /// A collection of cuts to apply on the physical phase space
  struct CutsList final : SteeredObject<CutsList> {
    explicit CutsList(const ParametersList&);

    static ParametersDescription description();

    void setParameters(const ParametersList&) override;

    friend std::ostream& operator<<(std::ostream&, const CutsList&);

    cuts::Initial initial;        ///< Cuts on the initial particles kinematics
    cuts::Central central;        ///< Cuts on the central system produced
    PerIdCuts central_particles;  ///< Cuts on the central individual particles
    cuts::Remnants remnants;      ///< Cuts on the beam remnants system
  };
}  // namespace cepgen

#endif
