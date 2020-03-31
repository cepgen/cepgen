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
      /// PDG ids of all known particles
      /// \note From \cite Beringer:1900zz :
      /// > The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.
      enum PdgId : pdgid_t
      {
        invalid = 0,
        down = 1, up = 2,
        electron = 11, muon = 13, tau = 15,
        gluon = 21, photon = 22, W = 24,
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

      /// Add a new particle definition to the library
      void define( const ParticleProperties& props );
      const std::vector<pdgid_t> particles() const; ///< All particles ids in this library
      void dump() const; ///< Dump all particles in this library
      size_t size() const; ///< Number of particles defined in this library

      //--- per-particles information

      bool has( pdgid_t ) const;
      /// All physical properties for one particle
      const ParticleProperties& operator()( pdgid_t ) const;
      const std::string& name( pdgid_t ) const; ///< Human-readable name for this particle
      double colours( pdgid_t ) const; ///< Colour factor for this particle
      double mass( pdgid_t ) const; ///< Particle mass (in GeV)
      double width( pdgid_t ) const; ///< Resonance width (in GeV)
      double charge( pdgid_t ) const; ///< Electric charge (in \f$e\f$) for this particle

    private:
      explicit PDG();
      /** \note Indexing variable: PDG id of particle */
      std::unordered_map<pdgid_t,ParticleProperties> particles_;
  };
}

#endif

