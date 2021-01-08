#include "CepGen/Physics/Cuts.h"

#include <iostream>
#include <algorithm>

namespace cepgen
{
  namespace cuts
  {
    //--------------------------------------------------------------------
    // physics system kinematic properties
    //--------------------------------------------------------------------

    Central::Central()
    {
      limits_.resize( CentralProperties::num_properties );
      limits_.at( e_pt_single ) = Property{ "pt", "Single particle pt (GeV/c)", Limits() };
      limits_.at( e_eta_single ) = Property{ "eta", "Single particle eta", Limits() };
      limits_.at( e_rapidity_single ) = Property{ "rapidity", "Single particle rapidity", Limits() };
      limits_.at( e_energy_single ) = Property{ "energy", "Single particle energy (GeV)", Limits() };
      limits_.at( e_mass_single ) = Property{ "mass", "Single particle mass (GeV/c^2)", Limits() };
      limits_.at( e_pt_sum ) = Property{ "ptsum", "System pt (GeV/c)", Limits() };
      limits_.at( e_eta_sum ) = Property{ "etasum", "System eta", Limits() };
      limits_.at( e_energy_sum ) = Property{ "energysum", "System energy (GeV)", Limits() };
      limits_.at( e_mass_sum ) = Property{ "invmass", "System mass (GeV/c^2)", Limits() };
      limits_.at( e_pt_diff ) = Property{ "ptdiff", "System D(pt) (GeV/c)", Limits() };
      limits_.at( e_phi_diff ) = Property{ "dphi", "System D(phi) (rad)", Limits() };
      limits_.at( e_rapidity_diff ) = Property{ "rapiditydiff", "System D(Y)", Limits() };
    }

    Initial::Initial()
    {
      limits_.resize( InitialProperties::num_properties );
      limits_.at( e_q2 ) = Property{ "q2", "Virtuality (GeV^2)", Limits() };
      limits_.at( e_qt ) = Property{ "qt", "Transverse virtuality (GeV)", Limits() };
      limits_.at( e_phi_qt ) = Property{ "phiqt", "Partons D(phi) (rad)", Limits() };
    }

    Remnants::Remnants()
    {
      limits_.resize( RemnantsProperties::num_properties );
      limits_.at( e_mx ) = Property{ "mx", "Diffractive mass (GeV/c^2)", Limits() };
      limits_.at( e_yj ) = Property{ "yj", "Diffractive jet rapidity", Limits() };
      limits_.at( e_xi ) = Property{ "xi", "Longit. fractional momentum loss", Limits() };
    }
  }

  //--------------------------------------------------------------------
  // utilitaries
  //--------------------------------------------------------------------

  std::vector<Cuts::Property>
  Cuts::list() const
  {
    std::vector<Property> out;
    std::copy_if( limits_.begin(), limits_.end(),
      std::back_inserter( out ),
      []( const Property& lim ) { return lim.limits.valid(); } );
    return out;
  }

  std::ostream&
  operator<<( std::ostream& os, const Cuts::Property& prop )
  {
    return os << "{" << prop.name << ": " << prop.limits << "}";
  }
}
