#include "CepGen/Physics/PDG.h"

namespace cepgen
{
  std::ostream&
  operator<<( std::ostream& os, const PDG& pc )
  {
    return os << std::string( PDGInfo::get()( pc ).human_name );
  }

  PDGInfo::PDGInfo()
  {
    add( PDG::invalid, { "[...]", 0, -1, -1., -1., 0., false } );
    //--- SM quarks
    add( PDG::down, { "down", "d", 3, 0.0048, 0., -1., true } );
    add( PDG::up, { "up", "u", 3, 0.0023, 0., 2., true } );
    add( PDG::strange, { "strange", "s", 3, 0.095, 0., -1., true } );
    add( PDG::charm, { "charm", "c", 3, 1.29, 0., 2., true } );
    add( PDG::bottom, { "bottom", "b", 3, 4.18, 0., -1., true } );
    add( PDG::top, { "top", "t", 3, 172.44, 0., 2., true } );
    //--- SM leptons
    add( PDG::electron, { "electron", "e± ", 1, 0.510998928e-3, 0., 3., true } );
    add( PDG::electronNeutrino, { "nu_e", "ν_e ", 1, 0., 0., 0., true } );
    add( PDG::muon, { "muon", "µ±  ", 1, 0.1056583715, 0., 3., true } );
    add( PDG::muonNeutrino, { "nu_mu", "ν_µ  ", 1, 0., 0., 0., true } );
    add( PDG::tau, { "tau", "τ±  ", 1, 1.77682, 0., 3., true } );
    add( PDG::muonNeutrino, { "nu_mu", "ν_τ  ", 1, 0., 0., 0., true } );
    //--- SM bosons
    add( PDG::gluon, { "gluon", "gluon", 9, 0., 0., 0., false } );
    add( PDG::photon, { "photon", "ɣ ", 0, 0., 0., 0., false } );
    add( PDG::Z, { "Z", "Z", 0, 91.1876, 2.4952, 0., false } );
    add( PDG::W, { "W", "W± ", 0, 80.385, 2.085, 3., false } );
    //--- nucleons
    add( PDG::proton, { "proton", "proton", 0, 0.938272046, 0., 3., false } );
    add( PDG::diffractiveProton, { "diff_proton", "diffr.proton", 0, 0., 0., 3., false } );
    add( PDG::neutron, { "neutron", "neutron", 0, 0.939565346, 0., 0., false } );
    //--- general mesons & baryons
    add( PDG::piPlus, { "pi_plus", "π±  ", 1, 0.13957018, -1., 3., false } );
    add( PDG::piZero, { "pi_zero", "π⁰  ", 1, 0.1349766, -1., 0., false } );
    add( PDG::KPlus, { "K_plus", "K± ", 1, 0.493677, -1., 3., false } );
    add( PDG::DPlus, { "D_plus", "D± ", 1, 1.86962, -1., 3., false } );
    add( PDG::rho770_0, { "rho770_0", "ρ(770)₀  ", 1, 0.77526, 0.150, 0., false } );
    add( PDG::rho1450_0, { "rho1450_0", "ρ(1450)₀  ", 1, 1.465, 0.400, 0., false } );
    add( PDG::rho1700_0, { "rho1700_0", "ρ(1700)₀  ", 1, 1.720, 0.250, 0., false } );
    add( PDG::h1380_1, { "h1380_1", "h(1380)₁ ", 1, 1.38619, false } );
    add( PDG::eta, { "eta", "η meson", 1, 0.547862, -1., 0., false } );
    add( PDG::omega782, { "omega782", "ω(782) ", 1, -1., 0., false } );
    add( PDG::Jpsi, { "Jpsi", "J/ψ ", 1, 3.0969, 92.9e-6 /* FIXME */, 0., false } );
    add( PDG::phi1680, { "phi1680", "ɸ(1680) ", 1, -1., -1., 0., false } );
    add( PDG::Upsilon1S, { "Upsilon1S", "Υ(1S) ", 1, 9.46030, 54.02e-6, 0., false } );
    add( PDG::Upsilon2S, { "Upsilon2S", "Υ(2S) ", 1, 10.02326, 31.98e-6, 0., false } );
    add( PDG::Upsilon3S, { "Upsilon3S", "Υ(3S) ", 1, 10.3552, 20.32e-6, 0., false } );
    add( PDG::pomeron, { "pomeron", "IP", 0, 0., 0., 0., false } );
    add( PDG::reggeon, { "reggeon", "IR", 0, 0., 0., 0., false } );
  }

  PDGInfo&
  PDGInfo::get()
  {
    static PDGInfo instance;
    return instance;
  }

  const ParticleProperties&
  PDGInfo::operator()( const PDG& id ) const
  {
    return particles_.at( id );
  }

  void
  PDGInfo::add( const PDG& id, const ParticleProperties& props )
  {
    particles_[id] = props;
  }
}
