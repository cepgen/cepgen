#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"

#include <iomanip>

namespace cepgen
{
  /*std::ostream&
  operator<<( std::ostream& os, pdgid_t pc )
  {
    return os << std::string( PDG::get()( pc ).description );
  }*/

  PDG::PDG()
  {
    define( invalid, { "[...]", 0, -1, -1., -1., 0, false } );
    //--- SM quarks
    define( down, { "down", "d", 3, 0.0048, 0., -1, true } );
    define( up, { "up", "u", 3, 0.0023, 0., 2, true } );
    define( strange, { "strange", "s", 3, 0.095, 0., -1, true } );
    define( charm, { "charm", "c", 3, 1.29, 0., 2, true } );
    define( bottom, { "bottom", "b", 3, 4.18, 0., -1, true } );
    define( top, { "top", "t", 3, 172.44, 0., 2, true } );
    //--- SM leptons
    define( electron, { "electron", "e± ", 1, 0.510998928e-3, 0., 3, true } );
    define( electronNeutrino, { "nu_e", "ν_e ", 1, 0., 0., 0, true } );
    define( muon, { "muon", "µ±  ", 1, 0.1056583715, 0., 3, true } );
    define( muonNeutrino, { "nu_mu", "ν_µ  ", 1, 0., 0., 0, true } );
    define( tau, { "tau", "τ±  ", 1, 1.77682, 0., 3, true } );
    define( muonNeutrino, { "nu_mu", "ν_τ  ", 1, 0., 0., 0, true } );
    //--- SM bosons
    define( gluon, { "gluon", "gluon", 9, 0., 0., 0, false } );
    define( photon, { "photon", "ɣ ", 0, 0., 0., 0, false } );
    define( Z, { "Z", "Z", 0, 91.1876, 2.4952, 0, false } );
    define( W, { "W", "W± ", 0, 80.385, 2.085, 3, false } );
    //--- nucleons
    define( proton, { "proton", "proton", 0, 0.938272046, 0., 3, false } );
    define( diffractiveProton, { "diff_proton", "diffr.proton", 0, 0., 0., 3, false } );
    define( neutron, { "neutron", "neutron", 0, 0.939565346, 0., 0, false } );
    //--- general mesons & baryons
    define( piPlus, { "pi_plus", "π±  ", 1, 0.13957018, -1., 3, false } );
    define( piZero, { "pi_zero", "π⁰  ", 1, 0.1349766, -1., 0, false } );
    define( KPlus, { "K_plus", "K± ", 1, 0.493677, -1., 3, false } );
    define( DPlus, { "D_plus", "D± ", 1, 1.86962, -1., 3, false } );
    define( rho770_0, { "rho770_0", "ρ(770)₀  ", 1, 0.77526, 0.150, 0, false } );
    define( rho1450_0, { "rho1450_0", "ρ(1450)₀  ", 1, 1.465, 0.400, 0, false } );
    define( rho1700_0, { "rho1700_0", "ρ(1700)₀  ", 1, 1.720, 0.250, 0, false } );
    define( h1380_1, { "h1380_1", "h(1380)₁ ", 1, 1.38619, 0, false } );
    define( eta, { "eta", "η meson", 1, 0.547862, -1., 0, false } );
    define( omega782, { "omega782", "ω(782) ", 1, -1., 0, false } );
    define( Jpsi, { "Jpsi", "J/ψ ", 1, 3.0969, 92.9e-6 /* FIXME */, 0, false } );
    define( phi1680, { "phi1680", "ɸ(1680) ", 1, -1., -1., 0, false } );
    define( Upsilon1S, { "Upsilon1S", "Υ(1S) ", 1, 9.46030, 54.02e-6, 0, false } );
    define( Upsilon2S, { "Upsilon2S", "Υ(2S) ", 1, 10.02326, 31.98e-6, 0, false } );
    define( Upsilon3S, { "Upsilon3S", "Υ(3S) ", 1, 10.3552, 20.32e-6, 0, false } );
    define( pomeron, { "pomeron", "IP", 0, 0., 0., 0, false } );
    define( reggeon, { "reggeon", "IR", 0, 0., 0., 0, false } );
  }

  PDG&
  PDG::get()
  {
    static PDG instance;
    return instance;
  }

  const ParticleProperties&
  PDG::operator()( int id ) const
  {
    return particles_.at( (pdgid_t)id );
  }

  void
  PDG::define( pdgid_t id, const ParticleProperties& props )
  {
    particles_[id] = props;
  }

  void
  PDG::dump() const
  {
    //--- first build a sorted vector out of the (unsorted) map
    std::vector<std::pair<int,ParticleProperties> > tmp;
    for ( const auto& prt : particles_ )
      tmp.emplace_back( prt.first, prt.second );
    std::sort( tmp.begin(), tmp.end(),
      []( const std::pair<int,ParticleProperties>& a, const std::pair<int,ParticleProperties>& b ) {
        return a.first < b.first;
      } );
    //--- then the proper dump begins
    std::ostringstream oss;
    for ( const auto& prt : tmp )
      if ( prt.first != PDG::invalid )
        oss
          << "\n  [" << std::setw( 8 ) << prt.first << "] "
          << std::setw( 12 ) << prt.second.name
          << " (" << prt.second.description << ")";
    CG_INFO( "PDG" ) << "List of particles registered:" << oss.str();
  }
}
