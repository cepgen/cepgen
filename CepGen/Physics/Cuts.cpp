#include "CepGen/Physics/Cuts.h"

#include <algorithm>
#include <iostream>

#include "CepGen/Core/ParametersList.h"

namespace cepgen {
  Cuts::Cuts(const ParametersList& params) {
    for (auto& cut : rawList())
      cut.limits = params.get<Limits>(cut.name);
  }

  namespace cuts {
    //--------------------------------------------------------------------
    // physics system kinematic properties
    //--------------------------------------------------------------------

    Central::Central() : Cuts(ParametersList()) { *this = Central(ParametersList()); }

    Central::Central(const ParametersList& params) : Cuts(params) {
      limits_.resize(CentralProperties::num_properties);
      limits_.at(e_pt_single) = Property{"pt", "Single particle pt (GeV/c)", params.get<Limits>("pt")};
      limits_.at(e_eta_single) = Property{"eta", "Single particle eta", params.get<Limits>("eta")};
      limits_.at(e_rapidity_single) = Property{"rapidity", "Single particle rapidity", params.get<Limits>("rapidity")};
      limits_.at(e_energy_single) = Property{"energy", "Single particle energy (GeV)", params.get<Limits>("energy")};
      limits_.at(e_mass_single) = Property{"mass", "Single particle mass (GeV/c^2)", params.get<Limits>("mass")};
      limits_.at(e_pt_sum) = Property{"ptsum", "System pt (GeV/c)", params.get<Limits>("ptsum")};
      limits_.at(e_eta_sum) = Property{"etasum", "System eta", params.get<Limits>("etasum")};
      limits_.at(e_energy_sum) = Property{"energysum", "System energy (GeV)", params.get<Limits>("energysum")};
      limits_.at(e_mass_sum) = Property{"invmass", "System mass (GeV/c^2)", params.get<Limits>("invmass")};
      limits_.at(e_pt_diff) = Property{"ptdiff", "System D(pt) (GeV/c)", params.get<Limits>("ptdiff")};
      limits_.at(e_phi_diff) = Property{"dphi", "System D(phi) (rad)", params.get<Limits>("dphi")};
      limits_.at(e_rapidity_diff) = Property{"rapiditydiff", "System D(Y)", params.get<Limits>("rapiditydiff")};
    }

    Initial::Initial(const ParametersList& params) : Cuts(params) {
      limits_.resize(InitialProperties::num_properties);
      limits_.at(e_q2) = Property{"q2", "Virtuality (GeV^2)", params.get<Limits>("q2")};
      limits_.at(e_qt) = Property{"qt", "Transverse virtuality (GeV)", params.get<Limits>("qt")};
      limits_.at(e_phi_qt) = Property{"phiqt", "Partons D(phi) (rad)", params.get<Limits>("phiqt")};
    }

    Remnants::Remnants(const ParametersList& params) : Cuts(params) {
      limits_.resize(RemnantsProperties::num_properties);
      limits_.at(e_mx) = Property{"mx", "Diffractive mass (GeV/c^2)", params.get<Limits>("mx")};
      limits_.at(e_yj) = Property{"yj", "Diffractive jet rapidity", params.get<Limits>("yj")};
      limits_.at(e_xi) = Property{"xi", "Longit. fractional momentum loss", params.get<Limits>("xi")};
    }
  }  // namespace cuts

  //--------------------------------------------------------------------
  // utilitaries
  //--------------------------------------------------------------------

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
