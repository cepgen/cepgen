#ifndef GenericCanvas_h
#define GenericCanvas_h

#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"

/**
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 25 Jul 2015
 */
class GenericCanvas : public TCanvas
{
  public:
    inline GenericCanvas() :
      TCanvas( "null" ), fBuilt( false ),
      fLegend( 0 ), fLegendX(.55), fLegendY( .74 ), fLegendNumEntries( 0 ),
      fUpperLabel( 0 ), fLabelsDrawn(false) {;}
    inline GenericCanvas( TString name, unsigned int width=500, unsigned int height=500, TString upper_label="" ) :
      TCanvas( name, "", width, height ), fBuilt( false ), fWidth( width ), fHeight( height ),
      fLegend( 0 ), fLegendX( .55 ), fLegendY( .74 ), fLegendNumEntries( 0 ),
      fUpperLabelText( upper_label ), fUpperLabel( 0 ), fLabelsDrawn( false ) { Build(); }
    inline GenericCanvas( TString name, TString upper_label ) :
      TCanvas( name, "", 500, 500 ), fBuilt( false ), fWidth( 500 ), fHeight( 500 ),
      fLegend( 0 ), fLegendX( .55 ), fLegendY( .74 ), fLegendNumEntries( 0 ),
      fUpperLabelText( upper_label ), fUpperLabel( 0 ), fLabelsDrawn( false ) { Build(); }

    inline virtual ~GenericCanvas() {
      if ( fLegend ) delete fLegend;
      if ( fUpperLabel ) delete fUpperLabel;
    }

    inline void SetUpperLabel( TString text="" ) {
      if ( text=="" ) return;
      fUpperLabelText = text;
      fUpperLabel = new TPaveText( .5, .945, .945, .98, "ndc" );
      fUpperLabel->SetMargin( 0. );
      fUpperLabel->SetFillColor( kWhite );
      fUpperLabel->SetLineColor( kWhite );
      fUpperLabel->SetLineWidth( 0 );
      fUpperLabel->SetShadowColor( kWhite );
      fUpperLabel->SetTextFont( 43 );
      fUpperLabel->SetTextAlign( 33 );
      fUpperLabel->SetTextSize( 18 );
      fUpperLabel->AddText( fUpperLabelText );
      fUpperLabel->Draw( "same" );
    }

    inline void Save( TString ext="png", TString path="." ) {
      bool valid_ext = true;
      valid_ext |= ( strcmp( ext.Data(), "png" )!=0 );
      valid_ext |= ( strcmp( ext.Data(), "pdf" )!=0 );
      if ( !valid_ext ) return;
      DrawLabels();
      //printf( "File saved as %s\n", Form( "%s/%s.%s", path.Data(), TCanvas::GetName(), ext.Data() ) );
      TCanvas::SaveAs( Form( "%s/%s.%s", path.Data(), TCanvas::GetName(), ext.Data() ) );
      //c1->SetLogz();
      //TCanvas::SaveAs( Form( "%s/%s_logscale.%s", path.Data(), TCanvas::GetName(), ext.Data() ) );
    }
    inline TPad* Pad() { return c1; }

    void AddLegendEntry( const TObject* obj_, TString label_, Option_t* option_="lf" ) {
      fLegend->AddEntry( obj_, label_, option_ );
      fLegendNumEntries++;
      if ( fLegendNumEntries>3 ) {
        fLegend->SetY1( fLegend->GetY1()-( fLegendNumEntries-3 )*0.015 );
      }
    }
   
    inline void Prettify( TH1* o ) {
      Prettify( o->GetXaxis(), o->GetYaxis() );
      o->SetTitle("");
    }
    inline void Prettify( TF1* o ) {
      Prettify( o->GetXaxis(), o->GetYaxis() );
      o->SetTitle("");
    }
    inline void Prettify( TGraph* o ) {
      Prettify( o->GetHistogram()->GetXaxis(), o->GetHistogram()->GetYaxis() );
      o->SetTitle("");
    }
    inline void Prettify( TAxis* x, TAxis* y=0 ) {
      // O so pretty axis!
      x->SetTitleFont( 43 );
      x->SetTitleSize( 28 );
      x->SetLabelFont( 43 );
      x->SetLabelSize( 22 );
      x->SetTitleOffset( 0.85 );
      
      if ( !y ) return;
      y->SetTitleFont( 43 );
      y->SetTitleSize( 28 );
      y->SetLabelFont( 43 );
      y->SetLabelSize( 22 );
      y->SetTitleOffset( 1.18 );
    }

  protected:
    inline void Build() {
      if ( fBuilt ) return;
      fLegend = new TLegend( fLegendX, fLegendY, fLegendX+.35, fLegendY+.15 );
      fLegend->SetFillColor( kWhite );
      fLegend->SetLineColor( kWhite );
      fLegend->SetLineWidth( 0 );
      fLegend->SetTextFont( 43 );
      fLegend->SetTextSize( 22 );
      DrawGrid();
      fBuilt = true;
    }
    
  private:
    inline void DrawLabels() {
      //if (fLabelsDrawn) return;
      
      if ( fLegend->GetNRows()!=0 ) fLegend->Draw();
      SetUpperLabel( fUpperLabelText );
      fLabelsDrawn = true;

      gStyle->SetMarkerStyle( 20 );
      gStyle->SetMarkerSize( .87 );
      gStyle->SetTitleFont( 43, "XYZ" );
      gStyle->SetTitleSize( 24, "XYZ" );
      //gStyle->SetTitleOffset( 2., "Y" );
      gStyle->SetLabelFont( 43, "XYZ" );
      gStyle->SetLabelSize( 20, "XY" );
      gStyle->SetLabelSize( 15, "Z" );
      gStyle->SetTitleOffset( 0.9, "X" );
      gStyle->SetTitleOffset( 1.1, "Y" );
      gStyle->SetHistLineColor( kBlack );
      gStyle->SetHistLineWidth( 2 );          
    }
    
    inline void DrawGrid() {
      TCanvas::cd();
      gStyle->SetOptStat( 0 );

      TCanvas::Divide( 1, 2 );
      c1 = dynamic_cast<TPad*>( TCanvas::GetPad( 1 ) );
      c2 = dynamic_cast<TPad*>( TCanvas::GetPad( 2 ) );
      c1->SetPad( 0., 0., 1., 1. );
      c2->SetPad( 0., 0., 1., 0. );
      c1->SetBottomMargin( 0.12 );
      c1->SetLeftMargin( 0.16 );
      c1->SetRightMargin( 0.05 );
      c1->SetTopMargin( 0.08 );
      TCanvas::cd( 1 );
     
      c1->SetTicks( 1, 1 );
      //c1->SetGrid( 1, 1 );
    }

    bool fBuilt;
    TPad *c1, *c2;
    double fWidth, fHeight;
    TLegend *fLegend;
    double fLegendX, fLegendY;
    unsigned int fLegendNumEntries;
    TPaveText *fLabel1, *fLabel2;
    TString fUpperLabelText;
    TPaveText *fUpperLabel;
    bool fLabelsDrawn;
};

#endif
