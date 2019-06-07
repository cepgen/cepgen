#ifndef CepGen_Physics_PDG_h
#define CepGen_Physics_PDG_h

#include <string>
#include <unordered_map>

namespace cepgen
{
  /** \brief PDG ids of all known particles
   * \note From \cite Beringer:1900zz :
   * `The Monte Carlo particle numbering scheme [...] is intended to facilitate interfacing between event generators, detector simulators, and analysis packages used in particle physics.`
   */
  enum class PDG
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
  std::ostream& operator<<( std::ostream& os, const PDG& pc );

  class PDGInfo
  {
    public:
      /// Retrieve a unique instance of this factory
      static PDGInfo& get() {
        static PDGInfo instance;
        return instance;
      }
      /// Default destructor
      ~PDGInfo() = default;
      struct ParticleProperties
      {
        std::string name, human_name;
        short colours;
        double mass, width, charge;
        bool isFermion;
      };

      void add( const PDG& id, const ParticleProperties& props );
      const ParticleProperties& operator()( const PDG& ) const;

    private:
      explicit PDGInfo();
      /** \note Indexing variable: PDG id of particle */
      std::unordered_map<PDG,ParticleProperties> particles_;
  };
}

#endif

