#ifndef CepGen_Physics_PDG_h
#define CepGen_Physics_PDG_h

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/ParticleProperties.h"

#include <unordered_map>

namespace cepgen
{
  /// A singleton holding all physics constants associated to particles
  class PDG
  {
    public:
      /** \brief PDG ids of all known particles
       * \note From \cite Beringer:1900zz :
       * `The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.`
       */
      enum PdgId : pdgid_t
      {
        invalid = 0,
        //--- fundamental particles
        down = 1, up = 2, strange = 3, charm = 4, bottom = 5, top = 6,
        electron = 11, electronNeutrino = 12,
        muon = 13, muonNeutrino = 14,
        tau = 15, tauNeutrino = 16,
        gluon = 21, photon = 22, Z = 23, W = 24,
        //--- composite particles
        piPlus = 211, piZero = 111,
        KPlus = 321, DPlus = 411,
        rho770_0 = 113, rho1450_0 = 100113, rho1700_0 = 30113,
        eta = 221, omega782 = 223,
        h1380_1 = 10333,
        Jpsi= 443,
        phi1680 = 100333,
        Upsilon1S = 553, Upsilon2S = 100553, Upsilon3S = 200553,
        proton = 2212, neutron = 2112,
        pomeron = 990, reggeon = 110,
        diffractiveProton = 9902210
      };

      /// Retrieve a unique instance of this particles info collection
      static PDG& get();
      PDG( const PDG& ) = delete;
      void operator=( const PDG& ) = delete;
      /// Default destructor
      ~PDG() = default;

      void define( pdgid_t id, const ParticleProperties& props );
      const ParticleProperties& operator()( pdgid_t ) const;
      void dump() const;
      std::string name( pdgid_t ) const;

    private:
      explicit PDG();
      /** \note Indexing variable: PDG id of particle */
      std::unordered_map<pdgid_t,ParticleProperties> particles_;
  };
}

#endif

