#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "THStack.h"
#include "../TreeEvent.h"
#include "../Canvas.h"
#include "CepGen/Physics/Particle.h"

void
validation()
{
  vector<pair<string, vector<string> > > files = {
    { "elastic", { "output_original_lpair_elastic_pt25.root", "output_cepgen_lpair_elastic_pt25.root" } },
    { "single diss.", { "output_original_lpair_singlediss_pt25.root", "output_cepgen_lpair_singlediss_pt25.root" } },
    { "double diss.", { "output_original_lpair_doublediss_pt25.root", "output_cepgen_lpair_doublediss_pt25.root" } },
  };

  map<string,TH1D*> m_plt[6];
  for ( unsigned short i = 0; i < 6; ++i ) {
    m_plt[i] = {
      { "invm", new TH1D( Form( "invm_%i", i ), "Dilepton invariant mass\\d#sigma/dm\\GeV", 20, 0., 800. ) },
      { "ptpair", new TH1D( Form( "ptpair_%i", i ), "Dilepton p_{T}\\d#sigma/dp_{T}\\GeV", 40, 0., 160. ) },
      { "singlept", new TH1D( Form( "singlept_%i", i ), "Single lepton p_{T}\\d#sigma/dp_{T}\\GeV", 25, 25., 150. ) },
      { "singleeta", new TH1D( Form( "singleeta_%i", i ), "Single lepton #eta\\d#sigma/d#eta\\?.2f", 15, -2.5, 2.5 ) },
    };
  }

  CepGen::TreeEvent ev;
  unsigned int n = 0;
  for ( const auto& kin : files ) {
    const string kinematics = kin.first;
    for ( const auto& file : kin.second ) {
      cout << "reading " << file << endl;
      auto tree = dynamic_cast<TTree*>( TFile::Open( ( string( "samples/" )+file ).c_str() )->Get( "h4444" ) );
      if ( !tree ) return;
      ev.attach( tree );
      const unsigned long long num_entries = tree->GetEntriesFast();
      double weight = 1./num_entries;
      for ( unsigned int i = 0; i < num_entries; ++i ) {
        tree->GetEntry( i );
        if ( i==0 ) weight *= ev.xsect;
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
        m_plt[n]["invm"]->Fill( ( lep1+lep2 ).M(), weight );
        m_plt[n]["ptpair"]->Fill( ( lep1+lep2 ).Pt(), weight );
        m_plt[n]["singlept"]->Fill( lep1.Pt(), weight );
        m_plt[n]["singleeta"]->Fill( lep1.Eta(), weight );
      }
      n++;
    }
  }

  //----- plotting part

  Color_t cols[] = { kBlack, kBlue+1, kRed+1 };

  for ( const auto& plt : m_plt[0] ) {
    CepGen::Canvas c( Form( "valid_%s", plt.first.c_str() ), "" );
    
    THStack hs;
    TH1D* h_sum = (TH1D*)m_plt[0][plt.first]->Clone();
    h_sum->Clear();
    for ( unsigned short i = 0; i < 6; ++i ) {
      TH1D* plot = m_plt[i][plt.first];
      h_sum->Add( plot );
      const char* style = ( i%2 == 0 ) ? "hist" : "e1";
      if ( i%2 == 0 ) {
        plot->SetLineWidth( 2 );
        plot->SetLineStyle( 1+int( i/2 ) );
        plot->SetLineColor( cols[i/2] );
      }
      else {
        plot->SetLineColor( kBlack );
        plot->SetLineWidth( 2 );
        plot->SetMarkerStyle( 20+int( i/3 ) );
        plot->SetMarkerColor( cols[i/2] );
        c.AddLegendEntry( plot, files[i/2].first.c_str() );
      }
      hs.Add( plot, style );
    }
    hs.Draw( "nostack" );
    hs.SetMaximum( h_sum->GetMaximum()/2. );
    hs.SetTitle( m_plt[0][plt.first]->GetTitle() );
    c.Prettify( hs.GetHistogram() );
    hs.SetTitle( "" );
    c.SetLogy();
    c.Save( "pdf" );
  }
}
