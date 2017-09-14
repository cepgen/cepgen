#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"

#include "../TreeEvent.h"
#include "../Canvas.h"

#include "CepGen/Physics/Particle.h"

void
simple_reader( const char* file="output.root" )
{
  map<string,TH1D*> m_plt = {
    { "invm", new TH1D( "invm", "Pair invariant mass\\d#sigma/dm\\GeV", 100, 0., 1000. ) },
    { "ptpair", new TH1D( "ptpair", "Pair p_{T}\\d#sigma/dp_{T}\\GeV?.1f", 100, 0., 1000. ) },
    { "singlept", new TH1D( "singlept", "Single p_{T}\\d#sigma/dp_{T}\\GeV", 125, 25., 150. ) },
    { "singleeta", new TH1D( "singleeta", "Single #eta\\d#sigma/d#eta\\?.2f", 5, -2.5, 2.5 ) },
    { "acopl", new TH1D( "acopl", "Acoplanarity |#Delta#phi/#pi|\\d#sigma/d#phi\\?.2f", 100, 0., 1. ) },
    { "mx", new TH1D( "mx", "Dissociated proton mass\\d#sigma/dM_{X}\\GeV", 100, 0., 1000. ) },
  };

  CepGen::TreeEvent ev;
  unsigned int n = 0;
  auto tree = dynamic_cast<TTree*>( TFile::Open( file )->Get( "h4444" ) );
  if ( !tree ) return;

  ev.attach( tree );
  const unsigned long long num_entries = tree->GetEntriesFast()/1;
  double weight = 1./num_entries;
  for ( unsigned int i = 0; i < num_entries; ++i ) {
    tree->GetEntry( i );
    if ( i==0 ) weight *= ev.xsect;
    TLorentzVector lep1, lep2, ip1, ip2, op1, op2;
    bool has_lepton1 = false;
    for ( unsigned short j = 0; j < ev.np ; ++j ) {
      switch ( ev.role[j] ) {
        case CepGen::Particle::IncomingBeam1: ip1.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::IncomingBeam2: ip2.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::OutgoingBeam1: op1.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::OutgoingBeam2: op2.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); break;
        case CepGen::Particle::CentralSystem: {
          if ( !has_lepton1 ) { lep1.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); has_lepton1 = true; }
          else                { lep2.SetPtEtaPhiM( ev.pt[j], ev.eta[j], ev.phi[j], ev.M[j] ); }
        } break;
        default: break;
      }
    }
    m_plt["invm"]->Fill( ( lep1+lep2 ).M(), weight );
    m_plt["ptpair"]->Fill( ( lep1+lep2 ).Pt(), weight );
    m_plt["singlept"]->Fill( lep1.Pt(), weight );
    m_plt["singleeta"]->Fill( lep1.Eta(), weight );
    m_plt["acopl"]->Fill( lep1.DeltaPhi( lep2 )/M_PI, weight );
    m_plt["mx"]->Fill( op1.M(), weight );
  }

  //----- plotting part

  for ( const auto& plt : m_plt ) {
    CepGen::Canvas c( Form( "cepgen_%s", plt.first.c_str() ), "CepGen simulation, pp at #sqrt{s} = 13 TeV" );

    TH1D* plot = m_plt[plt.first];
    plot->SetLineWidth( 2 );
    plot->Draw( "e1" );

    //CepGen::PaveText label( 0.75, 0.5, 0.8, 0.85, "My beautiful box" );
    //label.SetTextSize( 0.04 );
    //label.Draw();
    c.Prettify( plot );
    c.SetGrid();
    c.SetLogy();
    c.Save( "pdf" );
  }
}
