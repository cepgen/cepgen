#ifndef Test_TreeInfo_h
#define Test_TreeInfo_h

#include "TFile.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include <string>

namespace CepGen
{
  /// All useful information about a generation run
  struct TreeRun
  {
    double sqrt_s; ///< Centre of mass energy for beam particles
    double xsect; ///< Process cross section, in pb
    double errxsect; ///< Uncertainty on process cross section, in pb
    unsigned int num_events; ///< Number of events generated in run
    unsigned int litigious_events; ///< Number of litigious events in run
    /// ROOT tree used for storage/retrieval of this run information
    TTree* tree;

    TreeRun() : tree( NULL ) { clear(); }
    /// Reinitialise the run tree
    void clear() {
      sqrt_s = -1.;
      xsect = errxsect = -1.;
      num_events = litigious_events = 0;
    }
    /// Populate the run tree
    void create() {
      tree = new TTree( "run", "a tree containing information on the previous run" );
      if ( !tree ) return;
      tree->Branch( "xsect", &xsect, "xsect/D" );
      tree->Branch( "errxsect", &errxsect, "errxsect/D" );
      tree->Branch( "num_events", &num_events, "num_events/i" );
      tree->Branch( "litigious_events", &litigious_events, "litigious_events/i" );
      tree->Branch( "sqrt_s", &sqrt_s, "sqrt_s/D" );
    }
    /// Fill the run tree
    void fill() {
      tree->Fill();
    }
    /// Attach the run tree reader to a given file
    void attach( const char* filename, const char* run_tree = "run" ) {
      attach( TFile::Open( filename ), run_tree );
    }
    /// Attach the run tree reader to a given tree
    void attach( TFile* file, const char* run_tree = "run" ) {
      tree = dynamic_cast<TTree*>( file->Get( run_tree ) );
      if ( !tree ) return;
      tree->SetBranchAddress( "xsect", &xsect );
      tree->SetBranchAddress( "errxsect", &errxsect );
      tree->SetBranchAddress( "num_events", &num_events );
      tree->SetBranchAddress( "litigious_events", &litigious_events );
      tree->SetBranchAddress( "sqrt_s", &sqrt_s );
      if ( tree->GetEntriesFast() > 1 )
        std::cerr << "The run tree has more than one entry." << std::endl;
      tree->GetEntry( 0 );
    }
  };

  /// All useful information about a generated event
  struct TreeEvent
  {
    // book a sufficienly large number to allow the large multiplicity
    // of excited proton fragmentation products
    static const unsigned short maxpart = 5000; ///< Maximal number of particles in event
    /// Tree for which the event is booked
    TTree* tree;
    /// A pointer to the file opened for storage/retrieval
    std::unique_ptr<TFile> file;

    float gen_time; ///< Event generation time
    float tot_time; ///< Total event generation time
    int nremn_ch[2], nremn_nt[2];
    int np; ///< Number of particles in the event
    std::vector<ROOT::Math::XYZTVector> momentum, *pMom;
    double pt[maxpart]; ///< Particles transverse momentum
    double eta[maxpart]; ///< Particles pseudo-rapidity
    double phi[maxpart]; ///< Particles azimutal angle
    double rapidity[maxpart]; ///< Particles rapidity
    double E[maxpart]; ///< Particles energy, in GeV
    double m[maxpart]; ///< Particles mass, in GeV/c\f${}^2\f$
    double charge[maxpart]; ///< Particles charges, in e
    int pdg_id[maxpart]; ///< Integer particles PDG id
    int parent1[maxpart]; ///< First particles mother
    int parent2[maxpart]; ///< Last particles mother
    int stable[maxpart]; ///< Whether the particle must decay or not
    int role[maxpart]; ///< Particles role in the event
    int status[maxpart]; ///< Integer status code

    TreeEvent() : tree( nullptr ), pMom( nullptr ) {
      clear();
    }
    /// Reinitialise the event content
    void clear() {
      gen_time = tot_time = 0.;
      for ( unsigned short i = 0; i < 2; ++i )
        nremn_ch[i] = nremn_nt[i] = 0;
      np = 0;
      momentum.clear();
      for ( unsigned short i = 0; i < maxpart; ++i ) {
        pt[i] = eta[i] = phi[i] = rapidity[i] = E[i] = m[i] = charge[i] = 0.;
        pdg_id[i] = parent1[i] = parent2[i] = stable[i] = role[i] = status[i] = 0;
      }
    }
    /// Fill the tree with a new event
    void fill() {
      if ( !tree )
        throw std::runtime_error( "TreeEvent: Trying to fill a non-existent tree!" );

      tree->Fill();
      clear();
    }
    /// Populate the tree and all associated branches
    void create( TTree* t ) {
      tree = t;
      if ( !tree ) return;
      tree->Branch( "npart", &np, "npart/I" );
      tree->Branch( "nremn_charged", nremn_ch, "nremn_charged[2]/I" );
      tree->Branch( "nremn_neutral", nremn_nt, "nremn_neutral[2]/I" );
      tree->Branch( "role", role, "role[npart]/I" );
      pMom = &momentum;
      tree->Branch( "momentum", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &pMom );
      tree->Branch( "pt", pt, "pt[npart]/D" );
      tree->Branch( "eta", eta, "eta[npart]/D" );
      tree->Branch( "phi", phi, "phi[npart]/D" );
      tree->Branch( "rapidity", rapidity, "rapidity[npart]/D" );
      tree->Branch( "E", E, "E[npart]/D" );
      tree->Branch( "m", m, "m[npart]/D" );
      tree->Branch( "charge", charge, "charge[npart]/D" );
      tree->Branch( "pdg_id", pdg_id, "pdg_id[npart]/I" );
      tree->Branch( "parent1", parent1, "parent1[npart]/I" );
      tree->Branch( "parent2", parent2, "parent2[npart]/I" );
      tree->Branch( "stable", stable, "stable[npart]/I" );
      tree->Branch( "status", status, "status[npart]/I" );
      tree->Branch( "generation_time", &gen_time, "generation_time/F" );
      tree->Branch( "total_time", &tot_time, "total_time/F" );
    }
    /// Attach the event tree reader to a given file
    void attach( const char* filename, const char* events_tree = "events" ) {
      file.reset( TFile::Open( filename ) );
      attach( file.get(), events_tree );
    }
    /// Attach the event tree reader to a given ROOT file
    void attach( TFile* f, const char* events_tree = "events" ) {
      tree = dynamic_cast<TTree*>( f->Get( events_tree ) );
      attach( tree );
    }
    /// Attach the event tree reader to a given tree
    void attach( TTree* t ) {
      tree = t;
      if ( !tree ) return;
      tree->SetBranchAddress( "npart", &np );
      tree->SetBranchAddress( "nremn_charged", nremn_ch );
      tree->SetBranchAddress( "nremn_neutral", nremn_ch );
      tree->SetBranchAddress( "role", role );
      tree->SetBranchAddress( "momentum", &pMom );
      tree->SetBranchAddress( "pt", pt );
      tree->SetBranchAddress( "eta", eta );
      tree->SetBranchAddress( "phi", phi );
      tree->SetBranchAddress( "rapidity", rapidity );
      tree->SetBranchAddress( "E", E );
      tree->SetBranchAddress( "m", m );
      tree->SetBranchAddress( "charge", charge );
      tree->SetBranchAddress( "pdg_id", pdg_id );
      tree->SetBranchAddress( "parent1", parent1 );
      tree->SetBranchAddress( "parent2", parent2 );
      tree->SetBranchAddress( "stable", stable );
      tree->SetBranchAddress( "status", status );
      tree->SetBranchAddress( "generation_time", &gen_time );
      tree->SetBranchAddress( "total_time", &tot_time );
    }
  };
}

#endif


