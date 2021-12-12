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

#ifndef CepGen_Physics_Cuts_h
#define CepGen_Physics_Cuts_h

#include <unordered_map>
#include <vector>

#include "CepGen/Core/ParametersDescription.h"
#include "CepGen/Physics/ParticleProperties.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  class ParametersList;

  /// Constraints to be applied on the events kinematics
  class Cuts {
  public:
    /// Define a cut from parameters list
    explicit Cuts(const ParametersList&);

    static ParametersDescription description();

    /// Modify a few parameters values
    void setParameters(const ParametersList&);
    /// Retrieve the cuts list into a steering parameters list
    ParametersList parameters(bool full = false) const;

    /// A set of properties for a given cut
    struct Property {
      Property() = default;
      explicit Property(const std::string& name, const std::string& descr, const ParametersList&);
      std::string name, description;
      Limits limits;
    };
    friend std::ostream& operator<<(std::ostream&, const Property&);

    /// A collection of limits properties
    std::vector<Property>& rawList() { return limits_; }
    /// A collection of const-qualified limits properties
    std::vector<Property> list() const;

  protected:
    /// List of limits associated to this phase space cuts definition
    std::vector<Property> limits_;
  };

  /// A namespace for all kinematic cuts
  namespace cuts {
    /// Centrally produced particles phase space cuts
    class Central final : public Cuts {
    public:
      Central();
      explicit Central(const ParametersList&);

      /// single particle transverse momentum
      Limits& pt_single() { return limits_.at(e_pt_single).limits; }
      /// single particle transverse momentum
      const Limits& pt_single() const { return limits_.at(e_pt_single).limits; }
      /// single particle pseudo-rapidity
      Limits& eta_single() { return limits_.at(e_eta_single).limits; }
      /// single particle pseudo-rapidity
      const Limits& eta_single() const { return limits_.at(e_eta_single).limits; }
      /// single particle rapidity
      Limits& rapidity_single() { return limits_.at(e_rapidity_single).limits; }
      /// single particle rapidity
      const Limits& rapidity_single() const { return limits_.at(e_rapidity_single).limits; }
      /// single particle energy
      Limits& energy_single() { return limits_.at(e_energy_single).limits; }
      /// single particle energy
      const Limits& energy_single() const { return limits_.at(e_energy_single).limits; }
      /// single particle mass
      Limits& mass_single() { return limits_.at(e_mass_single).limits; }
      /// single particle mass
      const Limits& mass_single() const { return limits_.at(e_mass_single).limits; }
      /// multiparticle system transverse momentum
      Limits& pt_sum() { return limits_.at(e_pt_sum).limits; }
      /// multiparticle system transverse momentum
      const Limits& pt_sum() const { return limits_.at(e_pt_sum).limits; }
      /// multiparticle system pseudo-rapidity
      Limits& eta_sum() { return limits_.at(e_eta_sum).limits; }
      /// multiparticle system pseudo-rapidity
      const Limits& eta_sum() const { return limits_.at(e_eta_sum).limits; }
      /// multiparticle system energy
      Limits& energy_sum() { return limits_.at(e_energy_sum).limits; }
      /// multiparticle system energy
      const Limits& energy_sum() const { return limits_.at(e_energy_sum).limits; }
      /// multiparticle system invariant mass
      Limits& mass_sum() { return limits_.at(e_mass_sum).limits; }
      /// multiparticle system invariant mass
      const Limits& mass_sum() const { return limits_.at(e_mass_sum).limits; }
      /// transverse momentum balance between the central particles
      Limits& pt_diff() { return limits_.at(e_pt_diff).limits; }
      /// transverse momentum balance between the central particles
      const Limits& pt_diff() const { return limits_.at(e_pt_diff).limits; }
      /// azimuthal angles difference between the central particles
      Limits& phi_diff() { return limits_.at(e_phi_diff).limits; }
      /// azimuthal angles difference between the central particles
      const Limits& phi_diff() const { return limits_.at(e_phi_diff).limits; }
      /// rapidity balance between the central particles
      Limits& rapidity_diff() { return limits_.at(e_rapidity_diff).limits; }
      /// rapidity balance between the central particles
      const Limits& rapidity_diff() const { return limits_.at(e_rapidity_diff).limits; }

    private:
      enum CentralProperties {
        e_pt_single,
        e_eta_single,
        e_rapidity_single,
        e_energy_single,
        e_mass_single,
        e_pt_sum,
        e_eta_sum,
        e_energy_sum,
        e_mass_sum,
        e_pt_diff,
        e_phi_diff,
        e_rapidity_diff,
        num_properties
      };
    };

    /// Initial parton-like particles phase space cuts
    class Initial final : public Cuts {
    public:
      explicit Initial(const ParametersList&);

      /// parton virtuality
      Limits& q2() { return limits_.at(e_q2).limits; }
      /// parton virtuality
      const Limits& q2() const { return limits_.at(e_q2).limits; }
      /// parton transverse virtuality
      Limits& qt() { return limits_.at(e_qt).limits; }
      /// parton transverse virtuality
      const Limits& qt() const { return limits_.at(e_qt).limits; }
      /// parton azimuthal angle difference
      Limits& phi_qt() { return limits_.at(e_phi_qt).limits; }
      /// parton azimuthal angle difference
      const Limits& phi_qt() const { return limits_.at(e_phi_qt).limits; }

    private:
      enum InitialProperties { e_q2, e_qt, e_phi_qt, num_properties };
    };

    /// Outgoing beam remnant-like particles phase space cuts
    class Remnants final : public Cuts {
    public:
      explicit Remnants(const ParametersList&);

      /// diffractive mass
      Limits& mx() { return limits_.at(e_mx).limits; }
      /// diffractive mass
      const Limits& mx() const { return limits_.at(e_mx).limits; }
      /// diffractive jet rapidity
      Limits& yj() { return limits_.at(e_yj).limits; }
      /// diffractive jet rapidity
      const Limits& yj() const { return limits_.at(e_yj).limits; }
      /// longitudinal momentum fraction
      Limits& xi() { return limits_.at(e_xi).limits; }
      /// longitudinal momentum fraction
      const Limits& xi() const { return limits_.at(e_xi).limits; }

    private:
      enum RemnantsProperties { e_mx, e_yj, e_xi, num_properties };
    };
  }  // namespace cuts
  /// Collection of cuts to be applied on all particle with a given PDG id
  typedef std::unordered_map<pdgid_t, cuts::Central> PerIdCuts;
}  // namespace cepgen

#endif
