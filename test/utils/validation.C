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
    { "Elastic", { "output_original_lpair_elastic_pt25.root", "output_cepgen_lpair_elastic_pt25.root" } },
    { "Single-dissociative", { "output_original_lpair_singlediss_pt25.root", "output_cepgen_lpair_singlediss_pt25.root" } },
    { "Double-dissociative", { "output_original_lpair_doublediss_pt25.root", "output_cepgen_lpair_doublediss_pt25.root" } },
  };

  map<string,TH1D*> m_plt[6];
  for ( unsigned short i = 0; i < 6; ++i ) {
    m_plt[i] = {
      { "invm", new TH1D( Form( "invm_%i", i ), "Dilepton invariant mass\\d#sigma/dm\\GeV", 50, 0., 500. ) },
      { "ptpair", new TH1D( Form( "ptpair_%i", i ), "Dilepton p_{T}\\d#sigma/dp_{T}\\GeV?.1f", 64, 0., 160. ) },
      { "singlept", new TH1D( Form( "singlept_%i", i ), "Single lepton p_{T}\\d#sigma/dp_{T}\\GeV", 25, 25., 150. ) },
      { "singleeta", new TH1D( Form( "singleeta_%i", i ), "Single lepton #eta\\d#sigma/d#eta\\?.2f", 20, -2.5, 2.5 ) },
      { "acopl", new TH1D( Form( "acopl_%i", i ), "Dilepton |#Delta#phi/#pi|\\d#sigma/d#phi\\?.2f", 50, 0., 1. ) },
      { "mx", new TH1D( Form( "mx_%i", i ), "Dissociated proton mass\\d#sigma/dM_{X}\\GeV", 50, 0., 1000. ) },
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
      const unsigned long long num_entries = tree->GetEntriesFast()/1;
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
        m_plt[n]["acopl"]->Fill( lep1.DeltaPhi( lep2 )/M_PI, weight );
        m_plt[n]["mx"]->Fill( op1.M(), weight );
      }
      n++;
    }
  }

  //----- plotting part

  Color_t cols[] = { kRed+1, kBlue, kGreen+2 };

  for ( const auto& plt : m_plt[0] ) {
    THStack hs;
    TH1D* h_sum = (TH1D*)m_plt[0][plt.first]->Clone();
    h_sum->Clear();
    TH1D* h_lpair[3], *h_cepgen[3];
    for ( unsigned short i = 0; i < 3; ++i ) {
      h_lpair[i] = (TH1D*)h_sum->Clone();
      h_cepgen[i] = (TH1D*)h_sum->Clone();
    }
    vector<pair<TH1D*, string> > legends;
    for ( unsigned short i = 0; i < 6; ++i ) {
      TH1D* plot = m_plt[i][plt.first];
      h_sum->Add( plot );
      if ( i%2 == 0 ) {
        plot->SetLineWidth( 2 );
        plot->SetLineStyle( 1+int( i/2 ) );
        plot->SetLineColor( cols[i/2] );
        h_lpair[i/3]->Add( plot );
        hs.Add( plot, "hist" );
      }
      else {
        plot->SetLineColor( kBlack );
        plot->SetLineWidth( 2 );
        plot->SetMarkerStyle( 24+( i+1 )/3 );
        plot->SetMarkerColor( cols[i/2] );
        legends.emplace_back( plot, files[i/2].first.c_str() );
        h_cepgen[(i-1)/3]->Add( plot );
        hs.Add( plot, "e1" );
      }
    }
    { // comparison plot
      CepGen::Canvas c( Form( "valid_cepgen-vs-lpair_%s", plt.first.c_str() ), "CepGen + LPAIR simulations, pp at #sqrt{s} = 13 TeV" );
      hs.Draw( "nostack" );
      hs.SetMaximum( hs.GetHistogram()->GetMaximum()*2.5 );
      hs.GetHistogram()->SetTitle( m_plt[0][plt.first]->GetTitle() );
      const double size_x = 0.2, size_y = 0.2;
      //double pos_x = ( plt.first == "acopl" ) ? 0.15 : 
      double pos_x = 0.19, pos_y = 0.72;
      c.SetLegendX1( 0.5 );
      if ( plt.first == "acopl" ) { c.SetLegendY1( 0.17 ); c.SetLegendX1( 0.2 ); }
      if ( plt.first == "ptpair" || plt.first == "singlept" || plt.first == "invm" ) pos_y = 0.17;
      if ( plt.first == "invm" ) pos_x = 0.24;
      if ( plt.first == "singleeta" ) c.SetLegendY1( 0.17 );
      for ( const auto& leg : legends ) {
        c.AddLegendEntry( leg.first, leg.second.c_str(), "lp" );
      }

      CepGen::PaveText label( pos_x, pos_y, pos_x+size_x, pos_y+size_y, "LPAIR #gamma#gamma#rightarrow#mu^{+}#mu^{-}\\p_{T}(single #mu^{#pm}) > 25 GeV\\#||{#eta(single #mu^{#pm})} < 2.5\\M_{X} < 1 TeV" );
      label.SetTextSize( 0.04 );
      label.Draw();
      c.Prettify( hs.GetHistogram() );
      hs.GetHistogram()->SetTitle( "" );
      c.SetLogy();
      c.Save( "pdf" );
    }
    { // ratio plot
      CepGen::Canvas c( Form( "valid_%s_ratio", plt.first.c_str() ) );
      for ( unsigned short i = 0; i < 3; ++i ) {
        h_cepgen[i]->Divide( h_lpair[i] );
        h_cepgen[i]->Draw( ( i > 0 ) ? "same" : "" );
        h_cepgen[i]->SetMarkerStyle( 24+i );
        h_cepgen[i]->SetMarkerColor( cols[i] );
      }
      h_cepgen[0]->SetTitle( m_plt[0][plt.first]->GetTitle() );
      h_cepgen[0]->GetYaxis()->SetRangeUser( 0.2, 1.8 );
      c.Prettify( h_cepgen[0] );
      c.Save( "pdf" );
    }
  }
}
