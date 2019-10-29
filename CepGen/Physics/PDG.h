#ifndef CepGen_Physics_PDG_h
#define CepGen_Physics_PDG_h

#include "CepGen/Physics/ParticleProperties.h"

#include <unordered_map>
#include <vector>

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
        down = 1, up = 2,
        electron = 11, muon = 13, tau = 15,
        gluon = 21, photon = 22, W = 23,
        pomeron = 990, reggeon = 110, piZero = 111, piPlus = 211, eta = 221,
        phi1680 = 100333,
        proton = 2212, diffractiveProton = 9902210
      };

      /// Retrieve a unique instance of this particles info collection
      static PDG& get();
      PDG( const PDG& ) = delete;
      void operator=( const PDG& ) = delete;
      /// Default destructor
      ~PDG() = default;

      void define( const ParticleProperties& props );
      const ParticleProperties& operator()( pdgid_t ) const;
      const std::vector<pdgid_t> particles() const;
      void dump() const;
      size_t size() const;
      const std::string& name( pdgid_t ) const;
      short colours( pdgid_t ) const;
      double mass( pdgid_t ) const;
      double width( pdgid_t ) const;
      double charge( pdgid_t ) const;

    private:
      explicit PDG();
      /** \note Indexing variable: PDG id of particle */
      std::unordered_map<pdgid_t,ParticleProperties> particles_;
  };
}

#endif

