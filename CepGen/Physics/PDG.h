#ifndef CepGen_Physics_PDG_h
#define CepGen_Physics_PDG_h

#include <iostream>

namespace CepGen
{
  /** Unique identifier for a particle type. From \cite Beringer:1900zz :
   * `The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.`
   * \brief PDG ids of all known particles
   */
  enum class PDG
  {
    invalid = 0,
    //--- fundamental particles
    TopQuark = 6,
    Electron = 11, ElectronNeutrino = 12,
    Muon = 13, MuonNeutrino = 14,
    Tau = 15, TauNeutrino = 16,
    Gluon = 21, Photon = 22, Z = 23, W = 24,
    //--- composite particles
    PiPlus = 211, PiZero = 111,
    KPlus = 321, DPlus = 411,
    Rho770_0 = 113, Rho1450_0 = 100113, Rho1700_0 = 30113,
    Eta = 221, Omega782 = 223,
    h1380_1 = 10333,
    JPsi= 443,
    Phi1680 = 100333,
    Upsilon1S = 553, Upsilon2S = 100553, Upsilon3S = 200553,
    Proton = 2212, Neutron = 2112,
    Pomeron = 990, Reggeon = 110,
    DiffractiveProton = 9902210
  };
  std::ostream& operator<<( std::ostream& os, const PDG& pc );
}

#endif

