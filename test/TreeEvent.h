#ifndef Test_TreeEvent_h
#define Test_TreeEvent_h

namespace CepGen
{
  class TreeEvent
  {
    public:
      static constexpr unsigned short maxpart = 10;

      double xsect, errxsect;
      double mx_p1, mx_p2;
      int np;
      int nremn_ch[2], nremn_nt[2];
      double pt[maxpart], eta[maxpart], phi[maxpart], rapidity[maxpart];
      double E[maxpart], M[maxpart], charge[maxpart];
      int PID[maxpart], parentid[maxpart], isstable[maxpart], role[maxpart], status[maxpart];
      TLorentzVector kinematics[maxpart];
      float gen_time, tot_time;
      int hadr_trials, litigious_events;

    public:
      TreeEvent() :
        xsect( -1.0 ), errxsect( -1.0 ), mx_p1( 0. ), mx_p2( 0. ), np( 0 ),
        nremn_ch{ 0, 0 }, nremn_nt{ 0, 0 },
        gen_time( -1.0 ), tot_time( -1.0 ),
        hadr_trials( 0 ), litigious_events( 0 )
      {
        for ( unsigned short i = 0; i < maxpart; ++i ) {
          pt[i] = eta[i] = phi[i] = rapidity[i] = E[i] = M[i] = charge[i] = 0.0;
          PID[i] = parentid[i] = isstable[i] = role[i] = status[i] = 0;
        }
      }
      void create( TTree* tree ) {
        if ( !tree ) return;
        tree->Branch( "xsect", &xsect, "xsect/D" );
        tree->Branch( "errxsect", &errxsect, "errxsect/D" );
        tree->Branch( "MX1", &mx_p1, "MX1/D" );
        tree->Branch( "MX2", &mx_p2, "MX2/D" );
        tree->Branch( "npart", &np, "npart/I" );
        tree->Branch( "nremn_charged", nremn_ch, "nremn_charged[2]/I" );
        tree->Branch( "nremn_neutral", nremn_nt, "nremn_neutral[2]/I" );
        //tree->Branch( "kinematics", kinematics, "TLorentzVector[npart]" );
        tree->Branch( "pt", pt, "pt[npart]/D" );
        tree->Branch( "eta", eta, "eta[npart]/D" );
        tree->Branch( "phi", phi, "phi[npart]/D" );
        tree->Branch( "rapidity", rapidity, "rapidity[npart]/D" );
        tree->Branch( "E", E, "E[npart]/D" );
        tree->Branch( "m", M, "m[npart]/D" );
        tree->Branch( "charge", charge, "charge[npart]/D" );
        tree->Branch( "icode", PID, "icode[npart]/I" );
        tree->Branch( "parent", parentid, "parent[npart]/I" );
        tree->Branch( "stable", isstable, "stable[npart]/I" );
        tree->Branch( "role", role, "role[npart]/I" );
        tree->Branch( "status", status, "status[npart]/I" );
        tree->Branch( "generation_time", &gen_time, "generation_time/F" );
        tree->Branch( "total_time", &tot_time, "total_time/F" );
        tree->Branch( "hadronisation_trials", &hadr_trials, "hadronisation_trials/I" );
      }
      void attach( TTree* tree ) {
        tree->SetBranchAddress( "xsect", &xsect );
        tree->SetBranchAddress( "errxsect", &errxsect );
        tree->SetBranchAddress( "MX1", &mx_p1 );
        tree->SetBranchAddress( "MX2", &mx_p2 );
        tree->SetBranchAddress( "npart", &np );
        tree->SetBranchAddress( "nremn_charged", nremn_ch );
        tree->SetBranchAddress( "nremn_neutral", nremn_ch );
        tree->SetBranchAddress( "pt", pt );
        tree->SetBranchAddress( "eta", eta );
        tree->SetBranchAddress( "phi", phi );
        tree->SetBranchAddress( "rapidity", rapidity );
        tree->SetBranchAddress( "E", E );
        tree->SetBranchAddress( "m", M );
        tree->SetBranchAddress( "charge", charge );
        tree->SetBranchAddress( "icode", PID );
        tree->SetBranchAddress( "parent", parentid );
        tree->SetBranchAddress( "stable", isstable );
        tree->SetBranchAddress( "role", role );
        tree->SetBranchAddress( "status", status );
        tree->SetBranchAddress( "generation_time", &gen_time );
        tree->SetBranchAddress( "total_time", &tot_time );
        tree->SetBranchAddress( "hadronisation_trials", &hadr_trials );
      }
  };
}

#endif


