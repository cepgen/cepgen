/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#ifndef CepGen_EventInterfaces_ROOTTreeInfo_h
#define CepGen_EventInterfaces_ROOTTreeInfo_h

#include <TFile.h>
#include <TTree.h>

#include <exception>
#include <iostream>
#include <string>

#include "CepGen/Event/Event.h"

namespace ROOT {
  /// All useful information about a generation run
  class CepGenRun {
  public:
    static constexpr const char* TREE_NAME = "run";  ///< Output tree name

    static CepGenRun load(TFile*, const std::string& run_tree = TREE_NAME);
    static CepGenRun load(const std::string&, const std::string& run_tree = TREE_NAME);

    double sqrt_s{-1.};                ///< Centre of mass energy for beam particles
    double xsect{-1.};                 ///< Process cross section, in pb
    double errxsect{-1.};              ///< Uncertainty on process cross section, in pb
    unsigned int num_events{0};        ///< Number of events generated in run
    unsigned int litigious_events{0};  ///< Number of litigious events in run
    std::string process_name;          ///< Unique name of the process generated in this run
    std::string process_parameters;    ///< Serialised process parameters

    explicit CepGenRun();
    /// Reinitialise the run tree
    void clear();
    /// Populate the run tree
    void create();
    /// Retrieve the ROOT tree
    TTree* tree() { return tree_.get(); }
    /// Fill the run tree
    void fill();
    /// Attach the run tree reader to a given file
    void attach(const char* filename, const char* run_tree = TREE_NAME) { attach(TFile::Open(filename), run_tree); }
    /// Attach the run tree reader to a given tree
    void attach(TFile* file, const char* run_tree = TREE_NAME);

  private:
    /// ROOT tree used for storage/retrieval of this run information
    std::shared_ptr<TTree> tree_;
  };

  /// All useful information about a generated event
  class CepGenEvent {
  public:
    // book a sufficiently large number to allow the large multiplicity
    // of excited proton fragmentation products
    static constexpr size_t MAX_PART = 5000;            ///< Maximal number of particles in event
    static constexpr const char* TREE_NAME = "events";  ///< Output tree name

    static CepGenEvent load(TFile*, const std::string& events_tree = TREE_NAME);
    static CepGenEvent load(const std::string&, const std::string& events_tree = TREE_NAME);

    cepgen::Event::EventMetadata metadata;
    cepgen::Event* event{nullptr};
    float gen_time{-1.};        ///< Event generation time
    float tot_time{-1.};        ///< Total event generation time
    float weight{-1.};          ///< Event weight
    int np{0};                  ///< Number of particles in the event
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

    CepGenEvent();
    /// Reinitialise the event content
    void clear();
    /// Retrieve the ROOT tree
    TTree* tree() { return tree_.get(); }
    /// Populate the tree and all associated branches
    void create();
    /// Attach the event tree reader to a given file
    void attach(const char* filename, const char* events_tree = TREE_NAME) {
      file_.reset(TFile::Open(filename));
      attach(file_.get(), events_tree);
    }
    /// Attach the event tree reader to a given ROOT file
    void attach(TFile* f, const char* events_tree = TREE_NAME) {
      //--- special constructor to avoid the memory to be cleared at destruction time
      tree_ = std::shared_ptr<TTree>(dynamic_cast<TTree*>(f->Get(events_tree)), [=](TTree*) {});
      attach();
    }
    /// Attach the event tree reader to a given tree
    void attach() {
      if (!tree_)
        throw std::runtime_error("Failed to attach to the events TTree!");
      tree_->SetBranchAddress("npart", &np);
      tree_->SetBranchAddress("role", role);
      tree_->SetBranchAddress("pt", pt);
      tree_->SetBranchAddress("eta", eta);
      tree_->SetBranchAddress("phi", phi);
      tree_->SetBranchAddress("rapidity", rapidity);
      tree_->SetBranchAddress("E", E);
      tree_->SetBranchAddress("m", m);
      tree_->SetBranchAddress("charge", charge);
      tree_->SetBranchAddress("pdg_id", pdg_id);
      tree_->SetBranchAddress("parent1", parent1);
      tree_->SetBranchAddress("parent2", parent2);
      tree_->SetBranchAddress("stable", stable);
      tree_->SetBranchAddress("status", status);
      tree_->SetBranchAddress("weight", &weight);
      tree_->SetBranchAddress("generation_time", &gen_time);
      tree_->SetBranchAddress("total_time", &tot_time);
      tree_->SetBranchAddress("metadata", &metadata);
      tree_->SetBranchAddress("event", &event);
      tree_attached_ = true;
    }

    //--- direct cepgen::Event I/O helpers

    /// Fill the tree with a new event
    void fill(const cepgen::Event&, bool compress = false);
    /// Read the next event in the file
    bool next(cepgen::Event&);

  private:
    /// Tree for which the event is booked
    std::shared_ptr<TTree> tree_;
    std::unique_ptr<TFile> file_;
    bool tree_attached_{false};
    unsigned long long num_read_events_{0ull};
  };
}  // namespace ROOT

#endif
