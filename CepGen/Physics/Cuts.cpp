#include "CepGen/Physics/Cuts.h"

namespace cepgen
{
  std::vector<std::pair<std::string,Limits> >
  Cuts::list() const
  {
    std::vector<std::pair<std::string,Limits> > out;
    if ( pt_single.valid() )
      out.emplace_back( "Single particle pt (GeV/c)", pt_single );
    if ( eta_single.valid() )
      out.emplace_back( "Single particle eta", eta_single );
    if ( rapidity_single.valid() )
     out.emplace_back( "Single particle rapidity", rapidity_single );
    if ( energy_single.valid() )
     out.emplace_back( "Single particle energy (GeV)", energy_single );
    if ( mass_single.valid() )
     out.emplace_back( "Single particle mass (GeV/c²)", mass_single );
    if ( pt_sum.valid() )
     out.emplace_back( "System pt (GeV/c)", pt_sum );
    if ( eta_sum.valid() )
     out.emplace_back( "System eta", eta_sum );
    if ( energy_sum.valid() )
     out.emplace_back( "System energy", energy_sum );
    if ( mass_sum.valid() )
     out.emplace_back( "System mass", mass_sum );
    if ( pt_diff.valid() )
     out.emplace_back( "System Δpt (GeV/c)", pt_diff );
    if ( phi_pt_diff.valid() )
     out.emplace_back( "System Δɸ", phi_pt_diff );
    if ( rapidity_diff.valid() )
     out.emplace_back( "System ΔY", rapidity_diff );
    if ( q2.valid() )
     out.emplace_back( "Virtuality range (GeV²)", q2 );
    if ( qt.valid() )
     out.emplace_back( "Transverse virtuality range (GeV)", qt );
    if ( phi_qt.valid() )
     out.emplace_back( "Partons Δɸ range", phi_qt );
    return out;
  }
}

