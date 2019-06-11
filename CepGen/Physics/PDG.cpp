#include "CepGen/Physics/PDG.h"
#include "CepGen/Core/Exception.h"

#include <iomanip>

namespace cepgen
{
  PDG::PDG()
  {
    define( invalid, { "[...]", "", 0, -1, -1., 0, false } );
    //--- SM quarks
    define( down, { "down", "d", 3, 0.0048, 0., -1, true } );
    define( up, { "up", "u", 3, 0.0023, 0., 2, true } );
    define( strange, { "strange", "s", 3, 0.095, 0., -1, true } );
    define( charm, { "charm", "c", 3, 1.29, 0., 2, true } );
    define( bottom, { "bottom", "b", 3, 4.18, 0., -1, true } );
    define( top, { "top", "t", 3, 172.44, 0., 2, true } );
    //--- SM leptons
    define( electron, { "electron", "e\u00B1", 1, 0.510998928e-3, 0., 3, true } );
    define( electronNeutrino, { "nu_e", "\u03BD_e", 1, 0., 0., 0, true } );
    define( muon, { "muon", "\u03BC\u00B1", 1, 0.1056583715, 0., 3, true } );
    define( muonNeutrino, { "nu_mu", "\u03BD_\u03BC", 1, 0., 0., 0, true } );
    define( tau, { "tau", "\u03C4\u00B1", 1, 1.77682, 0., 3, true } );
    define( tauNeutrino, { "nu_tau", "\u03BD_\u03C4", 1, 0., 0., 0, true } );
    //--- SM bosons
    define( gluon, { "gluon", "g", 9, 0., 0., 0, false } );
    define( photon, { "photon", "\u03B3", 0, 0., 0., 0, false } );
    define( Z, { "Z", "Z", 0, 91.1876, 2.4952, 0, false } );
    define( W, { "W", "W\u00B1", 0, 80.385, 2.085, 3, false } );
    //--- nucleons
    define( proton, { "proton", "p", 0, 0.938272046, 0., 3, false } );
    define( diffractiveProton, { "diff_proton", "p\u002A", 0, 0., 0., 3, false } );
    define( neutron, { "neutron", "n", 0, 0.939565346, 0., 0, false } );
    //--- general mesons & baryons
    define( piPlus, { "pi_plus", "\u03C0\u00B1", 1, 0.13957018, -1., 3, false } );
    define( piZero, { "pi_zero", "\u03C0\u2070", 1, 0.1349766, -1., 0, false } );
    define( KPlus, { "K_plus", "K\u00B1", 1, 0.493677, -1., 3, false } );
    define( DPlus, { "D_plus", "D\u00B1", 1, 1.86962, -1., 3, false } );
    define( rho770_0, { "rho770_0", "\u03C1(770)\u2080", 1, 0.77526, 0.150, 0, false } );
    define( rho1450_0, { "rho1450_0", "\u03C1(1450)\u2080", 1, 1.465, 0.400, 0, false } );
    define( rho1700_0, { "rho1700_0", "\u03C1(1700)\u2080", 1, 1.720, 0.250, 0, false } );
    define( h1380_1, { "h1380_1", "h(1380)\u2081", 1, 1.38619, 0, false } );
    define( eta, { "eta", "\u03B7", 1, 0.547862, -1., 0, false } );
    define( omega782, { "omega782", "\u03C9(782)", 1, -1., 0, false } );
    define( Jpsi, { "Jpsi", "J/\u03C8", 1, 3.0969, 92.9e-6 /* FIXME */, 0, false } );
    define( phi1680, { "phi1680", "\u03A6(1680)", 1, -1., -1., 0, false } );
    define( Upsilon1S, { "Upsilon1S", "\u03A5(1S)", 1, 9.46030, 54.02e-6, 0, false } );
    define( Upsilon2S, { "Upsilon2S", "\u03A5(2S)", 1, 10.02326, 31.98e-6, 0, false } );
    define( Upsilon3S, { "Upsilon3S", "\u03A5(3S)", 1, 10.3552, 20.32e-6, 0, false } );
    define( pomeron, { "pomeron", "\u2119", 0, 0., 0., 0, false } );
    define( reggeon, { "reggeon", "\u211D", 0, 0., 0., 0, false } );
  }

  //--------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const ParticleProperties& prop )
  {
    return os << prop.name << "{"
      << "desc=" << prop.description
      << ",colours=" << prop.colours
      << ",mass=" << prop.mass
      << ",width=" << prop.width
      << ",charge=" << prop.charge
      << ( prop.fermion ? ",fermion" : "" )
      << "}";
  }

  //--------------------------------------------------------------------

  PDG&
  PDG::get()
  {
    static PDG instance;
    return instance;
  }

  const ParticleProperties&
  PDG::operator()( pdgid_t id ) const
  {
    try {
      return particles_.at( id );
    } catch ( const std::out_of_range& ) {
      throw CG_FATAL( "PDG" )
        << "Failed to retrieve particle properties for PDG id " << id << "!";
    }
  }

  void
  PDG::define( pdgid_t id, const ParticleProperties& props )
  {
    CG_DEBUG( "PDG:define" ) << "Adding a new particle with "
      << "PDG id=" << std::setw( 8 ) << id << ", " << props;
    particles_[id] = props;
  }

  std::string
  PDG::name( pdgid_t id ) const
  {
    return operator()( id ).description;
  }

  void
  PDG::dump() const
  {
    //--- first build a sorted vector out of the (unsorted) map
    std::vector<std::pair<pdgid_t,ParticleProperties> > tmp;
    for ( const auto& prt : particles_ )
      tmp.emplace_back( prt.first, prt.second );
    std::sort( tmp.begin(), tmp.end(),
      []( const std::pair<pdgid_t,ParticleProperties>& a,
          const std::pair<pdgid_t,ParticleProperties>& b ) {
        return a.first < b.first;
      } );
    //--- then the proper dump begins
    std::ostringstream oss;
    for ( const auto& prt : tmp )
      if ( prt.first != PDG::invalid )
        oss << "\n" << prt.second;
    CG_INFO( "PDG" ) << "List of particles registered:" << oss.str();
  }
}
