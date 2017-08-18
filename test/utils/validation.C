#include "../TreeEvent.h"
#include "../Canvas.h"
#include "CepGen/Physics/Particle.h"

void
validation( const string& cepgen_file="events.root", const string& orig_file="" )
{
  CepGen::TreeEvent ev;
  auto tree_val = dynamic_cast<TTree*>( TFile::Open( cepgen_file.c_str() )->Get( "h4444" ) );
  if ( !tree_val ) return;

  string plots[] = { "invm", "ptpair", "singlept", "singleeta" };
  const unsigned short num_plots = sizeof( plots ) / sizeof( plots[0] );

  map<string,TH1D*> m_plt[2];
  for ( unsigned short i = 0; i < 2; ++i ) {
    m_plt[i].insert( make_pair( "invm", new TH1D( Form( "invm_%i", i ), "Dilepton invariant mass\\Events\\GeV", 100, 0., 1000. ) ) );
    m_plt[i].insert( make_pair( "ptpair", new TH1D( Form( "ptpair_%i", i ), "Dilepton p_{T}\\Events\\GeV", 100, 0., 10. ) ) );
    m_plt[i].insert( make_pair( "singlept", new TH1D( Form( "singlept_%i", i ), "Single lepton p_{T}\\Events\\GeV", 100, 0., 100. ) ) );
    m_plt[i].insert( make_pair( "singleeta", new TH1D( Form( "singleeta_%i", i ), "Single lepton p_{T}\\Events\\GeV", 60, -3., 3. ) ) );
  }

  ev.attach( tree_val );
  for ( unsigned long long i = 0; i < tree_val->GetEntries(); ++i ) {
    tree_val->GetEntry( i );
    TLorentzVector lep1, lep2, ip1, ip2, op1, op2;
    for ( unsigned short j = 0; j < ev.np ; ++j ) {
      switch ( ev.role[j] ) {
        case CepGen::Particle::IncomingBeam1: ip1.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::IncomingBeam2: ip2.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::OutgoingBeam1: op1.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::OutgoingBeam2: op2.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::CentralParticle1: lep1.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::CentralParticle2: lep2.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        default: break;
      }
    }
    m_plt[0]["invm"]->Fill( ( lep1+lep2 ).M() );
    m_plt[0]["ptpair"]->Fill( ( lep1+lep2 ).Pt() );
    m_plt[0]["singlept"]->Fill( lep1.Pt() );
    m_plt[0]["singleeta"]->Fill( lep1.Eta() );
  }
  
  if ( orig_file != "" ) {
  }

  //----- plotting part

  for ( const auto& plt : plots ) {
    CepGen::Canvas c( Form( "valid_%s", plt.c_str() ), "" );
    THStack hs;
    hs.Add( m_plt[0][plt] );
    hs.Add( m_plt[1][plt] );
    hs.Draw( "nostack" );
    hs.SetTitle( m_plt[0][plt]->GetTitle() );
    c.Prettify( hs.GetHistogram() );
    hs.SetTitle( "" );
    c.Save( "pdf" );
  }
}
