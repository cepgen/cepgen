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

#include <algorithm>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Physics/Cuts.h"
#include "CepGen/Physics/Kinematics.h"

namespace cepgen {
  Cuts::Cuts(const ParametersList& params) { setParameters(params); }

  void Cuts::setParameters(const ParametersList& params) {
    for (auto& cut : rawList())
      params.fill<Limits>(cut.name, cut.limits);
  }

  ParametersList Cuts::parameters(bool full) const {
    ParametersList params;
    for (const auto& lim : list()) {
      params.set<Limits>(lim.name, lim.limits);
      if (full)
        params.set<double>(lim.name + "min", lim.limits.min()).set<double>(lim.name + "max", lim.limits.max());
    }
    return params;
  }

  std::vector<Cuts::Property> Cuts::list() const {
    std::vector<Property> out;
    // filtered list of valid limits/cuts
    std::copy_if(limits_.begin(), limits_.end(), std::back_inserter(out), [](const Property& lim) {
      return lim.limits.valid();
    });
    return out;
  }

  ParametersDescription Cuts::description() {
    auto desc = ParametersDescription();
    return desc;
  }

  namespace cuts {
    //--------------------------------------------------------------------
    // physics system kinematic properties
    //--------------------------------------------------------------------

    Central::Central() : Cuts(ParametersList()) { *this = Central(ParametersList()); }

    Central::Central(const ParametersList& params) : Cuts(params) {
      limits_.resize(CentralProperties::num_properties);
      limits_.at(e_pt_single) = Property("pt", "Single particle pt (GeV/c)", params);
      limits_.at(e_eta_single) = Property("eta", "Single particle eta", params);
      limits_.at(e_rapidity_single) = Property("rapidity", "Single particle rapidity", params);
      limits_.at(e_energy_single) = Property("energy", "Single particle energy (GeV)", params);
      limits_.at(e_mass_single) = Property("mass", "Single particle mass (GeV/c^2)", params);
      limits_.at(e_pt_sum) = Property("ptsum", "System pt (GeV/c)", params);
      limits_.at(e_eta_sum) = Property("etasum", "System eta", params);
      limits_.at(e_energy_sum) = Property("energysum", "System energy (GeV)", params);
      limits_.at(e_mass_sum) = Property("invmass", "System mass (GeV/c^2)", params);
      limits_.at(e_pt_diff) = Property("ptdiff", "System D(pt) (GeV/c)", params);
      limits_.at(e_phi_diff) = Property("dphi", "System D(phi) (rad)", params);
      limits_.at(e_rapidity_diff) = Property("rapiditydiff", "System D(Y)", params);
    }

    bool Central::contain(const Particles& parts, const Event*) const {
      Momentum mom_sum;
      for (const auto& part : parts) {
        const auto& mom = part.momentum();
        if (!pt_single().contains(mom.pt()))
          return false;
        if (!eta_single().contains(mom.eta()))
          return false;
        if (!rapidity_single().contains(mom.rapidity()))
          return false;
        if (!energy_single().contains(mom.energy()))
          return false;
        if (!mass_single().contains(mom.mass()))
          return false;
        mom_sum += mom;
      }
      if (!pt_sum().contains(mom_sum.pt()))
        return false;
      if (!eta_sum().contains(mom_sum.eta()))
        return false;
      if (!energy_sum().contains(mom_sum.energy()))
        return false;
      if (!mass_sum().contains(mom_sum.mass()))
        return false;
      if (parts.size() > 1) {  // look at correlations
        const auto &mom1 = parts.at(0).momentum(), &mom2 = parts.at(1).momentum();
        if (!pt_diff().contains(fabs(mom1.pt() - mom2.pt())))
          return false;
        if (!phi_diff().contains(mom1.deltaPhi(mom2)))
          return false;
        if (!rapidity_diff().contains(fabs(mom1.rapidity() - mom2.rapidity())))
          return false;
      }
      return true;
    }

    Initial::Initial(const ParametersList& params) : Cuts(params) {
      limits_.resize(InitialProperties::num_properties);
      limits_.at(e_q2) = Property("q2", "Virtuality (GeV^2)", params);
      limits_.at(e_qt) = Property("qt", "Transverse virtuality (GeV)", params);
      limits_.at(e_phi_qt) = Property("phiqt", "Partons D(phi) (rad)", params);
    }

    bool Initial::contain(const Particles&, const Event*) const {
      //if (!q2().contains(
      return true;
    }

    Remnants::Remnants(const ParametersList& params) : Cuts(params) {
      limits_.resize(RemnantsProperties::num_properties);
      limits_.at(e_mx) = Property("mx", "Diffractive mass (GeV/c^2)", params);
      limits_.at(e_yj) = Property("yj", "Diffractive jet rapidity", params);
      limits_.at(e_xi) = Property("xi", "Longit. fractional momentum loss", params);
    }

    bool Remnants::contain(const Particles& parts, const Event* evt) const {
      for (const auto& part : parts) {
        if (part.status() != Particle::Status::FinalState)
          return true;
        if (evt && !xi().contains(1. - part.momentum().pz() / (*evt)[*part.mothers().begin()].momentum().pz()))
          return false;
        if (!yj().contains(fabs(part.momentum().rapidity())))
          return false;
      }
      return true;
    }
  }  // namespace cuts

  //--------------------------------------------------------------------
  // utilitaries
  //--------------------------------------------------------------------

  Cuts::Property::Property(const std::string& name, const std::string& descr, const ParametersList& params)
      : name(name), description(descr) {
    if (params.has<Limits>(name))
      limits = params.get<Limits>(name);
  }

  std::ostream& operator<<(std::ostream& os, const Cuts::Property& prop) {
    return os << "{" << prop.name << ": " << prop.limits << "}";
  }

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  CutsList::CutsList(const ParametersList& params)
      : SteeredObject(params), initial(params_), central(params_), remnants(params_) {
    setParameters(params);
  }

  ParametersList CutsList::parameters(bool extended) const {
    ParametersList params;
    params += initial.parameters(extended);   // incoming partons
    params += central.parameters(extended);   // central particles
    params += remnants.parameters(extended);  // beam remnants
    // per-PDGid selection
    if (!central_particles.empty()) {
      ParametersList per_part;
      for (const auto& cuts_vs_part : central_particles)
        per_part.set<ParametersList>(std::to_string(cuts_vs_part.first), cuts_vs_part.second.parameters(extended));
      params.set<ParametersList>("cuts", per_part);
    }
    return params;
  }

  void CutsList::setParameters(const ParametersList& params) {
    SteeredObject::setParameters(params);
    initial.setParameters(params_);
    central.setParameters(params_);
    remnants.setParameters(params_);
    if (params_.has<Limits>("phiptdiff")) {
      CG_WARNING("Kinematics") << "\"phiptdiff\" parameter is deprecated! "
                               << "Please use \"phidiff\" instead.";
      params_.fill<Limits>("phiptdiff", central.phi_diff());  //legacy
    }
    if (params_.has<ParametersList>("cuts")) {  // per-particle cuts
      const auto& per_parts = steer<ParametersList>("cuts");
      for (const auto& part : per_parts.keys())
        central_particles[(pdgid_t)stoi(part)].setParameters(per_parts.get<ParametersList>(part));
    }

    //--- outgoing remnants
    // sanity check
    if (remnants.mx().min() < Kinematics::MX_MIN) {
      CG_WARNING("Kinematics:setParameters") << "Minimum diffractive mass range (" << remnants.mx() << ") is invalid. "
                                             << "It is now set to " << Kinematics::MX_MIN << " GeV/c^2.";
      remnants.mx().min() = Kinematics::MX_MIN;
    }
  }

  ParametersDescription CutsList::description() {
    auto desc = ParametersDescription();
    desc.add<Limits>("q2", {0., 1.e5}).setDescription("Incoming partons virtuality range");
    desc.add<double>("ptmin", 0.).setDescription("Outgoing central system object minimum transverse momentum (GeV/c)");
    desc.add<Limits>("mx", {Kinematics::MX_MIN, 1000.})
        .setDescription("Outgoing diffractive beam particles mass range (GeV/c^2)");
    return desc;
  }

  std::ostream& operator<<(std::ostream& os, const CutsList& cuts) {
    std::string sep;
    os << "initial: {";
    for (const auto& cut : cuts.initial.list())
      os << sep << cut, sep = ", ";
    os << "}, central: {";
    sep.clear();
    for (const auto& cut : cuts.central.list())
      os << sep << cut, sep = ", ";
    os << "}, remnants: {";
    sep.clear();
    for (const auto& cut : cuts.remnants.list())
      os << sep << cut, sep = ", ";
    return os << "}";
  }
}  // namespace cepgen
