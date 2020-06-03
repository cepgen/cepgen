#include "CepGen/Physics/Cuts.h"
#include <iostream>

namespace cepgen
{
  //--------------------------------------------------------------------
  // physics system kinematic properties
  //--------------------------------------------------------------------

  CentralCuts::CentralCuts()
  {
    limits_.resize( CentralCutsProperties::num_properties );
    limits_.at( e_pt_single ) = Property{ "pt", "Single particle pt (GeV/c)" };
    limits_.at( e_eta_single ) = Property{ "eta", "Single particle eta" };
    limits_.at( e_rapidity_single ) = Property{ "rapidity", "Single particle rapidity" };
    limits_.at( e_energy_single ) = Property{ "energy", "Single particle energy (GeV)" };
    limits_.at( e_mass_single ) = Property{ "mass", "Single particle mass (GeV/c^2)" };
    limits_.at( e_pt_sum ) = Property{ "ptsum", "System pt (GeV/c)" };
    limits_.at( e_eta_sum ) = Property{ "etasum", "System eta" };
    limits_.at( e_energy_sum ) = Property{ "energysum", "System energy (GeV)" };
    limits_.at( e_mass_sum ) = Property{ "invmass", "System mass (GeV/c^2)" };
    limits_.at( e_pt_diff ) = Property{ "ptdiff", "System D(pt) (GeV/c)" };
    limits_.at( e_phi_diff ) = Property{ "dphi", "System D(phi) (rad)" };
    limits_.at( e_rapidity_diff ) = Property{ "rapiditydiff", "System D(Y)" };
  }

  InitialCuts::InitialCuts()
  {
    limits_.resize( InitialCutsProperties::num_properties );
    limits_.at( e_q2 ) = Property{ "q2", "Virtuality (GeV^2)" };
    limits_.at( e_qt ) = Property{ "qt", "Transverse virtuality (GeV)" };
    limits_.at( e_phi_qt ) = Property{ "phiqt", "Partons D(phi) (rad)" };
  }

  RemnantsCuts::RemnantsCuts()
  {
    limits_.resize( RemnantsCutsProperties::num_properties );
    limits_.at( e_mx ) = Property{ "mx", "Diffractive mass (GeV/c^2)" };
    limits_.at( e_yj ) = Property{ "yj", "Diffractive jet rapidity" };
    limits_.at( e_xi ) = Property{ "xi", "Longit. fractional momentum loss" };
  }

  //--------------------------------------------------------------------
  // utilitaries
  //--------------------------------------------------------------------

  std::vector<Cuts::Property>
  Cuts::list()
  {
    std::vector<Property> out;
    for ( auto& lim : limits_ )
      if ( lim.limits.valid() )
        out.emplace_back( lim );
    return out;
  }

  std::vector<Cuts::Property>
  Cuts::list() const
  {
    std::vector<Property> out;
    for ( auto& lim : limits_ )
      if ( lim.limits.valid() )
        out.emplace_back( lim );
    return out;
  }

  std::ostream&
  operator<<( std::ostream& os, const Cuts::Property& prop )
  {
    return os << "{" << prop.name << ": " << prop.limits << "}";
  }
}
