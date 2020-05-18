#ifndef CepGen_Physics_Cuts_h
#define CepGen_Physics_Cuts_h

#include "CepGen/Physics/Limits.h"
#include "CepGen/Physics/ParticleProperties.h"

#include <vector>
#include <unordered_map>

namespace cepgen
{
  /// Constraints to be applied on the events kinematics
  class Cuts
  {
    public:
      explicit Cuts();

      /// A set of properties for a given cut
      struct Property
      {
        std::string name, description;
        Limits limits;
      };

      /// single particle transverse momentum
      Limits& pt_single() { return limits_.at( e_pt_single ).limits; }
      /// single particle transverse momentum
      const Limits& pt_single() const { return limits_.at( e_pt_single ).limits; }
      /// single particle pseudo-rapidity
      Limits& eta_single() { return limits_.at( e_eta_single ).limits; }
      /// single particle pseudo-rapidity
      const Limits& eta_single() const { return limits_.at( e_eta_single ).limits; }
      /// single particle rapidity
      Limits& rapidity_single() { return limits_.at( e_rapidity_single ).limits; }
      /// single particle rapidity
      const Limits& rapidity_single() const { return limits_.at( e_rapidity_single ).limits; }
      /// single particle energy
      Limits& energy_single() { return limits_.at( e_energy_single ).limits; }
      /// single particle energy
      const Limits& energy_single() const { return limits_.at( e_energy_single ).limits; }
      /// single particle mass
      Limits& mass_single() { return limits_.at( e_mass_single ).limits; }
      /// single particle mass
      const Limits& mass_single() const { return limits_.at( e_mass_single ).limits; }
      /// multiparticle system transverse momentum
      Limits& pt_sum() { return limits_.at( e_pt_sum ).limits; }
      /// multiparticle system transverse momentum
      const Limits& pt_sum() const { return limits_.at( e_pt_sum ).limits; }
      /// multiparticle system pseudo-rapidity
      Limits& eta_sum() { return limits_.at( e_eta_sum ).limits; }
      /// multiparticle system pseudo-rapidity
      const Limits& eta_sum() const { return limits_.at( e_eta_sum ).limits; }
      /// multiparticle system energy
      Limits& energy_sum() { return limits_.at( e_energy_sum ).limits; }
      /// multiparticle system energy
      const Limits& energy_sum() const { return limits_.at( e_energy_sum ).limits; }
      /// multiparticle system invariant mass
      Limits& mass_sum() { return limits_.at( e_mass_sum ).limits; }
      /// multiparticle system invariant mass
      const Limits& mass_sum() const { return limits_.at( e_mass_sum ).limits; }
      /// transverse momentum balance between the central particles
      Limits& pt_diff() { return limits_.at( e_pt_diff ).limits; }
      /// transverse momentum balance between the central particles
      const Limits& pt_diff() const { return limits_.at( e_pt_diff ).limits; }
      /// azimuthal angles difference between the central particles
      Limits& phi_diff() { return limits_.at( e_phi_diff ).limits; }
      /// azimuthal angles difference between the central particles
      const Limits& phi_diff() const { return limits_.at( e_phi_diff ).limits; }
      /// rapidity balance between the central particles
      Limits& rapidity_diff() { return limits_.at( e_rapidity_diff ).limits; }
      /// rapidity balance between the central particles
      const Limits& rapidity_diff()  const{ return limits_.at( e_rapidity_diff ).limits; }
      /// parton virtuality
      Limits& q2() { return limits_.at( e_q2 ).limits; }
      /// parton virtuality
      const Limits& q2() const { return limits_.at( e_q2 ).limits; }
      /// parton transverse virtuality
      Limits& qt() { return limits_.at( e_qt ).limits; }
      /// parton transverse virtuality
      const Limits& qt() const { return limits_.at( e_qt ).limits; }
      /// parton azimuthal angle difference
      Limits& phi_qt() { return limits_.at( e_phi_qt ).limits; }
      /// parton azimuthal angle difference
      const Limits& phi_qt() const { return limits_.at( e_phi_qt ).limits; }

      /// A collection of limits properties
      std::vector<Property>& rawList() { return limits_; }
      /// A collection of valid limits properties
      std::vector<Property> list();
      /// A collection of const-qualified limits properties
      std::vector<Property> list() const;

    private:
      enum Properties
      {
        e_pt_single, e_eta_single, e_rapidity_single, e_energy_single, e_mass_single,
        e_pt_sum, e_eta_sum, e_energy_sum, e_mass_sum,
        e_pt_diff, e_phi_diff, e_rapidity_diff,
        e_q2, e_qt, e_phi_qt,
        num_properties
      };
      std::vector<Property> limits_;
  };
  /// Collection of cuts to be applied on all particle with a given PDG id
  typedef std::unordered_map<pdgid_t,Cuts> PerIdCuts;
}

#endif
