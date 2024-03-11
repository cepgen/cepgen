#ifndef CepGenAddOns_PAWWrapper_PAWTreeInfo_h
#define CepGenAddOns_PAWWrapper_PAWTreeInfo_h

namespace paw
{
  /// All useful information about a generation run
  struct CepGenRun
  {
    static constexpr const char* TREE_NAME = "run"; ///< Output tree name
    float sqrt_s; ///< Centre of mass energy for beam particles
    float xsect; ///< Process cross section, in pb
    float errxsect; ///< Uncertainty on process cross section, in pb
    int num_events; ///< Number of events generated in run
    int litigious_events; ///< Number of litigious events in run

    CepGenRun() {
      clear();
    }
    /// Reinitialise the run tree
    void clear() {
      sqrt_s = -1.;
      xsect = errxsect = -1.;
      num_events = litigious_events = 0;
    }
  };

  /// All useful information about a generated event
  struct CepGenEvent
  {
    // book a sufficienly large number to allow the large multiplicity
    // of excited proton fragmentation products
    static constexpr size_t MAX_PART = 5000; ///< Maximal number of particles in event
    static constexpr const char* TREE_NAME = "events"; ///< Output tree name

    float gen_time; ///< Event generation time
    float tot_time; ///< Total event generation time
    float weight; ///< Event weight
    int np; ///< Number of particles in the event
    float pt[MAX_PART]; ///< Particles transverse momentum
    float eta[MAX_PART]; ///< Particles pseudo-rapidity
    float phi[MAX_PART]; ///< Particles azimutal angle
    float rapidity[MAX_PART]; ///< Particles rapidity
    float E[MAX_PART]; ///< Particles energy, in GeV
    float m[MAX_PART]; ///< Particles mass, in GeV/c\f${}^2\f$
    float charge[MAX_PART]; ///< Particles charges, in e
    int pdg_id[MAX_PART]; ///< Integer particles PDG id
    int parent1[MAX_PART]; ///< First particles mother
    int parent2[MAX_PART]; ///< Last particles mother
    int stable[MAX_PART]; ///< Whether the particle must decay or not
    int role[MAX_PART]; ///< Particles role in the event
    int status[MAX_PART]; ///< Integer status code

    CepGenEvent() {
      clear();
    }
    /// Reinitialise the event content
    void clear() {
      gen_time = tot_time = 0.;
      np = 0;
      for ( size_t i = 0; i < MAX_PART; ++i ) {
        pt[i] = eta[i] = phi[i] = rapidity[i] = E[i] = m[i] = charge[i] = 0.;
        pdg_id[i] = parent1[i] = parent2[i] = stable[i] = role[i] = status[i] = 0;
      }
    }
  };
}

#endif

