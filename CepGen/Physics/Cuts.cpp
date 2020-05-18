#include "CepGen/Physics/Cuts.h"

namespace cepgen
{
  Cuts::Cuts() :
    limits_( Properties::num_properties )
  {
    limits_.at( e_pt_single ) = Property{ "pt", "Single particle pt (GeV/c)" };
    limits_.at( e_eta_single ) = Property{ "eta", "Single particle eta" };
    limits_.at( e_rapidity_single ) = Property{ "rapidity", "Single particle rapidity" };
    limits_.at( e_energy_single ) = Property{ "energy", "Single particle energy (GeV)" };
    limits_.at( e_mass_single ) = Property{ "mass", "Single particle mass (GeV/c²)" };
    limits_.at( e_pt_sum ) = Property{ "ptsum", "System pt (GeV/c)" };
    limits_.at( e_eta_sum ) = Property{ "etasum", "System eta" };
    limits_.at( e_energy_sum ) = Property{ "energysum", "System energy (GeV)" };
    limits_.at( e_mass_sum ) = Property{ "invmass", "System mass (GeV/c²)" };
    limits_.at( e_pt_diff ) = Property{ "ptdiff", "System Δpt (GeV/c)" };
    limits_.at( e_phi_diff ) = Property{ "dphi", "System Δɸ (rad)" };
    limits_.at( e_rapidity_diff ) = Property{ "rapiditydiff", "System ΔY" };
    limits_.at( e_q2 ) = Property{ "q2", "Virtuality (GeV²)" };
    limits_.at( e_qt ) = Property{ "qt", "Transverse virtuality (GeV)" };
    limits_.at( e_phi_qt ) = Property{ "phiqt", "Partons Δɸ (rad)" };
  }

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
}

