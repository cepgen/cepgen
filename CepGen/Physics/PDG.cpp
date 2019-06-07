#include "CepGen/Physics/PDG.h"

namespace cepgen
{
  std::ostream&
  operator<<( std::ostream& os, const PDG& pc )
  {
    return os << PDGInfo::get()( pc ).human_name;
  }

  PDGInfo::PDGInfo()
  {
    add( PDG::invalid, ParticleProperties{ "[...]", 0, -1, -1., -1., 0., false } );
    //--- SM quarks
    add( PDG::down, ParticleProperties{ "down", "d", 3, 0.0048, 0., -1., true } );
    add( PDG::up, ParticleProperties{ "up", "u", 3, 0.0023, 0., 2., true } );
    add( PDG::strange, ParticleProperties{ "strange", "s", 3, 0.095, 0., -1., true } );
    add( PDG::charm, ParticleProperties{ "charm", "c", 3, 1.29, 0., 2., true } );
    add( PDG::bottom, ParticleProperties{ "bottom", "b", 3, 4.18, 0., -1., true } );
    add( PDG::top, ParticleProperties{ "top", "t", 3, 172.44, 0., 2., true } );
    //--- SM leptons
    add( PDG::electron, ParticleProperties{ "electron", "e± ", 1, 0.510998928e-3, 0., 3., true } );
    add( PDG::electronNeutrino, ParticleProperties{ "nu_e", "ν_e ", 1, 0., 0., 0., true } );
    add( PDG::muon, ParticleProperties{ "muon", "µ±  ", 1, 0.1056583715, 0., 3., true } );
    add( PDG::muonNeutrino, ParticleProperties{ "nu_mu", "ν_µ  ", 1, 0., 0., 0., true } );
    add( PDG::tau, ParticleProperties{ "tau", "τ±  ", 1, 1.77682, 0., 3., true } );
    add( PDG::muonNeutrino, ParticleProperties{ "nu_mu", "ν_τ  ", 1, 0., 0., 0., true } );
    //--- SM bosons
    add( PDG::gluon, ParticleProperties{ "gluon", "gluon", 9, 0., 0., 0., false } );
    add( PDG::photon, ParticleProperties{ "photon", "ɣ ", 0, 0., 0., 0., false } );
    add( PDG::Z, ParticleProperties{ "Z", "Z", 0, 91.1876, 2.4952, 0., false } );
    add( PDG::W, ParticleProperties{ "W", "W± ", 0, 80.385, 2.085, 3., false } );
    //--- mesons & baryons
    add( PDG::piPlus, ParticleProperties{ "pi_plus", "π±  ", 1, 0.13957018, -1., 3., false } );
    add( PDG::piZero, ParticleProperties{ "pi_zero", "π⁰  ", 1, 0.1349766, -1., 0., false } );
    add( PDG::KPlus, ParticleProperties{ "K_plus", "K± ", 1, 0.493677, -1., 3., false } );
    add( PDG::DPlus, ParticleProperties{ "D_plus", "D± ", 1, 1.86962, -1., 3., false } );
    add( PDG::rho770_0, ParticleProperties{ "rho770_0", "ρ(770)₀  ", 1, 0.77526, 0.150, 0., false } );
    add( PDG::rho1450_0, ParticleProperties{ "rho1450_0", "ρ(1450)₀  ", 1, 1.465, 0.400, 0., false } );
    add( PDG::rho1700_0, ParticleProperties{ "rho1700_0", "ρ(1700)₀  ", 1, 1.720, 0.250, 0., false } );
    add( PDG::h1380_1, ParticleProperties{ "h1380_1", "h(1380)₁ ", 1, 1.38619, false } );
    add( PDG::eta, ParticleProperties{ "eta", "η meson", 1, 0.547862, -1., 0., false } );
    add( PDG::omega782, ParticleProperties{ "omega782", "ω(782) ", 1, -1., -1., false } );
    add( PDG::Jpsi, ParticleProperties{ "Jpsi", "J/ψ ", 1, 3.0969, 92.9e-6 /* FIXME */, 0., false } );
    add( PDG::phi1680, ParticleProperties{ "phi1680", "ɸ(1680) ", 1, -1., -1., 0., false } );
    add( PDG::Upsilon1S, ParticleProperties{ "Upsilon1S", "Υ(1S) ", 1, 9.46030, 54.02e-6, 0., false } );
    add( PDG::Upsilon2S, ParticleProperties{ "Upsilon2S", "Υ(2S) ", 1, 10.02326, 31.98e-6, 0., false } );
    add( PDG::Upsilon3S, ParticleProperties{ "Upsilon3S", "Υ(3S) ", 1, 10.3552, 20.32e-6, 0., false } );
    add( PDG::proton, ParticleProperties{ "proton", "proton", 0, 0.938272046, 0., 3., false } );
    add( PDG::diffractiveProton, ParticleProperties{ "diff_proton", "diffr.proton", 0, 0., 0., 3., false } );
    add( PDG::neutron, ParticleProperties{ "neutron", "neutron", 0, 0.939565346, 0., 0., false } );
    add( PDG::pomeron, ParticleProperties{ "pomeron", "IP", 0, 0., 0., 0., false } );
    add( PDG::reggeon, ParticleProperties{ "reggeon", "IR", 0, 0., 0., 0., false } );
  }

  const PDGInfo::ParticleProperties&
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
