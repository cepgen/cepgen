#include "CepGen/Physics/Cuts.h"

#include <algorithm>
#include <iostream>

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  Cuts::Cuts(const ParametersList& params) { setParameters(params); }

  void Cuts::setParameters(const ParametersList& params) {
    for (auto& cut : rawList())
      params.fill<Limits>(cut.name, cut.limits);
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

    Initial::Initial(const ParametersList& params) : Cuts(params) {
      limits_.resize(InitialProperties::num_properties);
      limits_.at(e_q2) = Property("q2", "Virtuality (GeV^2)", params);
      limits_.at(e_qt) = Property("qt", "Transverse virtuality (GeV)", params);
      limits_.at(e_phi_qt) = Property("phiqt", "Partons D(phi) (rad)", params);
    }

    Remnants::Remnants(const ParametersList& params) : Cuts(params) {
      limits_.resize(RemnantsProperties::num_properties);
      limits_.at(e_mx) = Property("mx", "Diffractive mass (GeV/c^2)", params);
      limits_.at(e_yj) = Property("yj", "Diffractive jet rapidity", params);
      limits_.at(e_xi) = Property("xi", "Longit. fractional momentum loss", params);
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

  std::vector<Cuts::Property> Cuts::list() const {
    std::vector<Property> out;
    std::copy_if(limits_.begin(), limits_.end(), std::back_inserter(out), [](const Property& lim) {
      return lim.limits.valid();
    });
    return out;
  }

  std::ostream& operator<<(std::ostream& os, const Cuts::Property& prop) {
    return os << "{" << prop.name << ": " << prop.limits << "}";
  }
}  // namespace cepgen
