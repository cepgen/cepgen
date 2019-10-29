#ifndef CepGen_EventInterfaces_ROOTTreeInfo_h
#define CepGen_EventInterfaces_ROOTTreeInfo_h

#include "TFile.h"
#include "TTree.h"

#include <exception>
#include <string>
#include <iostream>

namespace cepgen { class Event; }

namespace ROOT
{
  /// All useful information about a generation run
  class CepGenRun
  {
    public:
      static constexpr const char* TREE_NAME = "run"; ///< Output tree name
      double sqrt_s; ///< Centre of mass energy for beam particles
      double xsect; ///< Process cross section, in pb
      double errxsect; ///< Uncertainty on process cross section, in pb
      unsigned int num_events; ///< Number of events generated in run
      unsigned int litigious_events; ///< Number of litigious events in run

      CepGenRun() {
        clear();
      }
      /// Reinitialise the run tree
      void clear() {
        sqrt_s = -1.;
        xsect = errxsect = -1.;
        num_events = litigious_events = 0;
      }
      /// Populate the run tree
      void create() {
        tree_ = std::make_shared<TTree>( TREE_NAME, "a tree containing information on the previous run" );
        if ( !tree_ )
          throw std::runtime_error( "Failed to create the run TTree!" );
        tree_->Branch( "xsect", &xsect, "xsect/D" );
        tree_->Branch( "errxsect", &errxsect, "errxsect/D" );
        tree_->Branch( "num_events", &num_events, "num_events/i" );
        tree_->Branch( "litigious_events", &litigious_events, "litigious_events/i" );
        tree_->Branch( "sqrt_s", &sqrt_s, "sqrt_s/D" );
      }
      /// Retrieve the ROOT tree
      TTree* tree() {
        return tree_.get();
      }
      /// Fill the run tree
      void fill() {
        tree_->Fill();
      }
      /// Attach the run tree reader to a given file
      void attach( const char* filename, const char* run_tree = TREE_NAME ) {
        attach( TFile::Open( filename ), run_tree );
      }
      /// Attach the run tree reader to a given tree
      void attach( TFile* file, const char* run_tree = TREE_NAME ) {
        //--- special constructor to avoid the memory to be cleared at destruction time
        tree_ = std::shared_ptr<TTree>( dynamic_cast<TTree*>( file->Get( run_tree ) ), [=]( TTree* ){} );
        if ( !tree_ )
          throw std::runtime_error( "Failed to attach to the run TTree!" );
        tree_->SetBranchAddress( "xsect", &xsect );
        tree_->SetBranchAddress( "errxsect", &errxsect );
        tree_->SetBranchAddress( "num_events", &num_events );
        tree_->SetBranchAddress( "litigious_events", &litigious_events );
        tree_->SetBranchAddress( "sqrt_s", &sqrt_s );
        if ( tree_->GetEntriesFast() > 1 )
          std::cerr << "The run tree has more than one entry." << std::endl;
        tree_->GetEntry( 0 );
      }

    private:
      /// ROOT tree used for storage/retrieval of this run information
      std::shared_ptr<TTree> tree_;
  };

  /// All useful information about a generated event
  class CepGenEvent
  {
    public:
      // book a sufficienly large number to allow the large multiplicity
      // of excited proton fragmentation products
      static constexpr size_t MAX_PART = 5000; ///< Maximal number of particles in event
      static constexpr const char* TREE_NAME = "events"; ///< Output tree name

      float gen_time; ///< Event generation time
      float tot_time; ///< Total event generation time
      float weight; ///< Event weight
      int nremn_ch[2], nremn_nt[2];
      int np; ///< Number of particles in the event
      double pt[MAX_PART]; ///< Particles transverse momentum
      double eta[MAX_PART]; ///< Particles pseudo-rapidity
      double phi[MAX_PART]; ///< Particles azimutal angle
      double rapidity[MAX_PART]; ///< Particles rapidity
      double E[MAX_PART]; ///< Particles energy, in GeV
      double m[MAX_PART]; ///< Particles mass, in GeV/c\f${}^2\f$
      double charge[MAX_PART]; ///< Particles charges, in e
      int pdg_id[MAX_PART]; ///< Integer particles PDG id
      int parent1[MAX_PART]; ///< First particles mother
      int parent2[MAX_PART]; ///< Last particles mother
      int stable[MAX_PART]; ///< Whether the particle must decay or not
      int role[MAX_PART]; ///< Particles role in the event
      int status[MAX_PART]; ///< Integer status code

      CepGenEvent() : tree_attached_( false ), num_read_events_( 0ull ) {
        clear();
      }
      /// Reinitialise the event content
      void clear() {
        gen_time = tot_time = 0.;
        for ( unsigned short i = 0; i < 2; ++i )
          nremn_ch[i] = nremn_nt[i] = 0;
        np = 0;
        for ( size_t i = 0; i < MAX_PART; ++i ) {
          pt[i] = eta[i] = phi[i] = rapidity[i] = E[i] = m[i] = charge[i] = 0.;
          pdg_id[i] = parent1[i] = parent2[i] = stable[i] = role[i] = status[i] = 0;
        }
      }
      /// Retrieve the ROOT tree
      TTree* tree() {
        return tree_.get();
      }
      /// Fill the tree with a new event
      void fill() {
        if ( !tree_ )
          throw std::runtime_error( "CepGenEvent: Trying to fill a non-existent tree!" );

        tree_->Fill();
        clear();
      }
      /// Populate the tree and all associated branches
      void create() {
        tree_ = std::make_shared<TTree>( TREE_NAME, "a tree containing information on events generated in previous run" );
        if ( !tree_ )
          throw std::runtime_error( "Failed to create the events TTree!" );
        tree_->Branch( "npart", &np, "npart/I" );
        tree_->Branch( "nremn_charged", nremn_ch, "nremn_charged[2]/I" );
        tree_->Branch( "nremn_neutral", nremn_nt, "nremn_neutral[2]/I" );
        tree_->Branch( "role", role, "role[npart]/I" );
        tree_->Branch( "pt", pt, "pt[npart]/D" );
        tree_->Branch( "eta", eta, "eta[npart]/D" );
        tree_->Branch( "phi", phi, "phi[npart]/D" );
        tree_->Branch( "rapidity", rapidity, "rapidity[npart]/D" );
        tree_->Branch( "E", E, "E[npart]/D" );
        tree_->Branch( "m", m, "m[npart]/D" );
        tree_->Branch( "charge", charge, "charge[npart]/D" );
        tree_->Branch( "pdg_id", pdg_id, "pdg_id[npart]/I" );
        tree_->Branch( "parent1", parent1, "parent1[npart]/I" );
        tree_->Branch( "parent2", parent2, "parent2[npart]/I" );
        tree_->Branch( "stable", stable, "stable[npart]/I" );
        tree_->Branch( "status", status, "status[npart]/I" );
        tree_->Branch( "weight", &weight, "weight/F" );
        tree_->Branch( "generation_time", &gen_time, "generation_time/F" );
        tree_->Branch( "total_time", &tot_time, "total_time/F" );
      }
      /// Attach the event tree reader to a given file
      void attach( const char* filename, const char* events_tree = TREE_NAME ) {
        file_.reset( TFile::Open( filename ) );
        attach( file_.get(), events_tree );
      }
      /// Attach the event tree reader to a given ROOT file
      void attach( TFile* f, const char* events_tree = TREE_NAME ) {
        //--- special constructor to avoid the memory to be cleared at destruction time
        tree_ = std::shared_ptr<TTree>( dynamic_cast<TTree*>( f->Get( events_tree ) ), [=]( TTree* ){} );
        attach();
      }
      /// Attach the event tree reader to a given tree
      void attach() {
        if ( !tree_ )
          throw std::runtime_error( "Failed to attach to the events TTree!" );
        tree_->SetBranchAddress( "npart", &np );
        tree_->SetBranchAddress( "nremn_charged", nremn_ch );
        tree_->SetBranchAddress( "nremn_neutral", nremn_ch );
        tree_->SetBranchAddress( "role", role );
        tree_->SetBranchAddress( "pt", pt );
        tree_->SetBranchAddress( "eta", eta );
        tree_->SetBranchAddress( "phi", phi );
        tree_->SetBranchAddress( "rapidity", rapidity );
        tree_->SetBranchAddress( "E", E );
        tree_->SetBranchAddress( "m", m );
        tree_->SetBranchAddress( "charge", charge );
        tree_->SetBranchAddress( "pdg_id", pdg_id );
        tree_->SetBranchAddress( "parent1", parent1 );
        tree_->SetBranchAddress( "parent2", parent2 );
        tree_->SetBranchAddress( "stable", stable );
        tree_->SetBranchAddress( "status", status );
        tree_->SetBranchAddress( "weight", &weight );
        tree_->SetBranchAddress( "generation_time", &gen_time );
        tree_->SetBranchAddress( "total_time", &tot_time );
        tree_attached_ = true;
      }

      //--- direct cepgen::Event I/O helpers

      /// Fill the tree with a new event
      void fill( const cepgen::Event&, bool compress = false );
      /// Read the next event in the file
      bool next( cepgen::Event& );

    private:
      /// Tree for which the event is booked
      std::shared_ptr<TTree> tree_;
      std::unique_ptr<TFile> file_;
      bool tree_attached_;
      unsigned long long num_read_events_;
  };
}

#endif

