/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGenRoot_ROOTTreeInfo_h
#define CepGenRoot_ROOTTreeInfo_h

#include <TFile.h>
#include <TTree.h>

#include "CepGen/Event/Event.h"

namespace ROOT {
  /// All useful information about a generation run
  class CepGenRun {
  public:
    static constexpr const char* TREE_NAME = "run";  ///< Output tree name

    static CepGenRun load(TFile*, const std::string& run_tree = TREE_NAME);
    static CepGenRun load(const std::string&, const std::string& run_tree = TREE_NAME);

    double sqrt_s{-1.};                ///< Centre of mass energy for beam particles
    double xsect{-1.};                 ///< Process cross-section, in pb
    double errxsect{-1.};              ///< Uncertainty on process cross-section, in pb
    unsigned int num_events{0};        ///< Events multiplicity generated in run
    unsigned int litigious_events{0};  ///< Litigious events multiplicity in run
    std::string process_name;          ///< Unique name of the process generated in this run
    std::string process_parameters;    ///< Serialised process parameters

    explicit CepGenRun();
    void clear();                                                       ///< Reinitialise the run tree
    void create();                                                      ///< Populate the run tree
    inline TTree* tree() const { return tree_.get(); }                  ///< Retrieve the ROOT tree
    void fill() const;                                                  ///< Fill the run tree
    void attach(TFile* file, const std::string& run_tree = TREE_NAME);  ///< Attach the run tree reader to a given tree
    /// Attach the run tree reader to a given file
    void attach(const std::string& filename, const std::string& run_tree = TREE_NAME) {
      attach(TFile::Open(filename.data()), run_tree);
    }

  private:
    std::shared_ptr<TTree> tree_;  ///< ROOT tree used for storage/retrieval of this run information
  };

  /// All useful information about a generated event
  class CepGenEvent {
  public:
    CepGenEvent() { clear(); }

    // book a large enough number to allow the large multiplicity of excited proton fragmentation products
    static constexpr size_t MAX_PART = 5000;            ///< Maximal particles multiplicity in event
    static constexpr const char* TREE_NAME = "events";  ///< Output tree name

    static CepGenEvent load(TFile*, const std::string& events_tree = TREE_NAME);
    static CepGenEvent load(const std::string&, const std::string& events_tree = TREE_NAME);

    cepgen::Event::EventMetadata metadata;
    float gen_time{-1.};        ///< Event generation time
    float tot_time{-1.};        ///< Total event generation time
    float weight{-1.};          ///< Event weight
    int np{0};                  ///< Particles multiplicity in the event
    double pt[MAX_PART];        ///< Particles transverse momentum
    double eta[MAX_PART];       ///< Particles pseudo-rapidity
    double phi[MAX_PART];       ///< Particles azimuthal angle
    double rapidity[MAX_PART];  ///< Particles rapidity
    double E[MAX_PART];         ///< Particles energy, in GeV
    double m[MAX_PART];         ///< Particles mass, in GeV/c\f${}^2\f$
    double charge[MAX_PART];    ///< Particles charges, in e
    int pdg_id[MAX_PART];       ///< Integer particles PDG id
    int parent1[MAX_PART];      ///< First particles mother
    int parent2[MAX_PART];      ///< Last particles mother
    int stable[MAX_PART];       ///< Whether the particle must decay or not
    int role[MAX_PART];         ///< Particles role in the event
    int status[MAX_PART];       ///< Integer status code

    void clear();                                ///< Reinitialise the event content
    TTree* tree() const { return tree_.get(); }  ///< Retrieve the ROOT tree
    void create();                               ///< Populate the tree and all associated branches

    void attach();                                                      ///< Attach the event tree reader to a tree
    void attach(TFile* f, const std::string& events_tree = TREE_NAME);  ///< Attach the event tree reader to a file
    /// Attach the event tree reader to a file
    void attach(const std::string& filename, const std::string& events_tree = TREE_NAME);

    // direct cepgen::Event I/O helpers

    void fill(const cepgen::Event&, bool compress = false);  ///< Fill the tree with a new event
    bool next(cepgen::Event&);                               ///< Read the next event in the file

  private:
    std::shared_ptr<TTree> tree_;  ///< Tree for which the event is booked
    std::unique_ptr<TFile> file_;
    bool tree_attached_{false};
    unsigned long long num_read_events_{0ull};
  };
}  // namespace ROOT

#endif
