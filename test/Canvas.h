#ifndef Canvas_h
#define Canvas_h

#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TStyle.h"

#include <string.h>

#define font_type(x) 130 + x

namespace cepgen {
  /// A "prettified" text box object
  class PaveText : public TPaveText {
  public:
    inline PaveText(const float& x1, const float& y1, const float& x2, const float& y2, const char* text = "")
        : TPaveText(x1, y1, x2, y2, "NDC") {
      TPaveText::SetTextAlign(13);
      if (strcmp(text, "") != 0) {
        TString txt = text;
        if (txt.Contains("\\")) {
          TObjArray* tok = txt.Tokenize("\\");
          for (int i = 0; i < tok->GetEntries(); ++i)
            TPaveText::AddText(dynamic_cast<TObjString*>(tok->At(i))->String());
        } else
          TPaveText::AddText(text);
      }
      TPaveText::SetFillColor(0);
      TPaveText::SetFillStyle(0);
      TPaveText::SetLineColor(0);
      TPaveText::SetLineWidth(0);
      TPaveText::SetShadowColor(0);
      TPaveText::SetTextFont(font_type(2));
      TPaveText::SetTextSize(0.058);
    }
  };

  /// A "prettified" generic figure canvas
  class Canvas : public TCanvas {
  public:
    /// Build a canvas from its name, title, and attributes
    /// \param[in] name Canvas name (and subsequently filename on save)
    /// \param[in] ratio Divide the canvas into a main and ratio plots subparts?
    explicit inline Canvas(const char* name, const char* title = "", bool ratio = false)
        :  //TCanvas( name, "", 450, 450 ),
          TCanvas(name, "", 600, 600),
          fTitle(title),
          fTopLabel(0),
          fLeg(0),
          fLegX1(0.5),
          fLegY1(0.75),
          fRatio(ratio) {
      gStyle->SetOptStat(0);
      Build();
    }
    inline ~Canvas() {
      if (fLeg)
        delete fLeg;
      if (fTopLabel)
        delete fTopLabel;
    }

    inline void SetSize(const float& size = 600) { TCanvas::SetCanvasSize(size, 600); }

    inline void Prettify(TH1* obj) {
      TAxis *x = dynamic_cast<TAxis*>(obj->GetXaxis()), *y = dynamic_cast<TAxis*>(obj->GetYaxis()),
            *z = dynamic_cast<TAxis*>(obj->GetZaxis());
      x->CenterTitle();
      y->CenterTitle();
      z->CenterTitle();
      x->SetLabelFont(font_type(3));
      x->SetLabelSize(20);
      x->SetTitleFont(font_type(3));
      x->SetTitleSize(29);
      y->SetLabelFont(font_type(3));
      y->SetLabelSize(20);
      y->SetTitleFont(font_type(3));
      y->SetTitleSize(29);
      z->SetLabelFont(font_type(3));
      z->SetLabelSize(16);
      z->SetTitleFont(font_type(3));
      z->SetTitleSize(29);
      if (fRatio) {
        x->SetTitleOffset(3.);
        x->SetLabelOffset(0.02);
      }
      y->SetTitleOffset(1.3);
      x->SetTickLength(0.03);
      y->SetTickLength(0.03);

      // axis titles
      TString ttle = obj->GetTitle();
      if (ttle.Contains("\\")) {
        TObjArray* tok = ttle.Tokenize("\\");
        TString x_title = "", y_title = "", unit = "", form_spec = "", distrib = "";
        if (tok->GetEntries() > 0)
          x_title = dynamic_cast<TObjString*>(tok->At(0))->String();
        if (tok->GetEntries() > 1)
          y_title = dynamic_cast<TObjString*>(tok->At(1))->String();
        if (tok->GetEntries() > 2) {
          unit = ((TObjString*)tok->At(2))->String();
          if (unit.Contains("?")) {  // extract format specifier
            TObjArray* tok2 = unit.Tokenize("?");
            if (tok2->GetEntries() > 1) {
              unit = dynamic_cast<TObjString*>(tok2->At(0))->String();
              form_spec = dynamic_cast<TObjString*>(tok2->At(1))->String();
            } else {
              unit = "";
              form_spec = dynamic_cast<TObjString*>(tok2->At(0))->String();
            }
          }
        }
        if (tok->GetEntries() > 3) {
          distrib = ((TObjString*)tok->At(3))->String();
        }
        if (!unit.IsNull() or !form_spec.IsNull()) {
          if (!unit.IsNull())
            x_title = Form("%s (%s)", x_title.Data(), unit.Data());
          if (!distrib.IsNull()) {
            if (!form_spec.IsNull()) {
              TString format = Form("%%s (%s / %%%s %%s)", distrib.Data(), form_spec.Data());
              y_title = Form(format.Data(), y_title.Data(), GetBinning(obj), unit.Data());
            } else
              y_title = Form("%s (%s / %d %s)",
                             y_title.Data(),
                             distrib.Data(),
                             static_cast<unsigned int>(GetBinning(obj)),
                             unit.Data());
          } else {
            if (!form_spec.IsNull()) {
              TString format = Form("%%s / %%%s %%s", form_spec.Data());
              y_title = Form(format.Data(), y_title.Data(), GetBinning(obj), unit.Data());
            } else
              y_title = Form("%s / %d %s", y_title.Data(), static_cast<unsigned int>(GetBinning(obj)), unit.Data());
          }
        }
        obj->GetXaxis()->SetTitle(x_title);
        obj->GetYaxis()->SetTitle(y_title);
        obj->SetTitle("");
      }
      //else obj->GetXaxis()->SetTitle(ttle);
    }

    inline void DrawDiagonal(const TH1* obj) {
      TLine l;
      l.SetLineWidth(2);
      l.SetLineColor(kGray);
      l.SetLineStyle(2);
      l.DrawLine(obj->GetXaxis()->GetXmin(),
                 obj->GetYaxis()->GetXmin(),
                 obj->GetXaxis()->GetXmax(),
                 obj->GetYaxis()->GetXmax());
    }

    inline void RatioPlot(TH1* obj1, const TH1* obj2, const TH1* obj3, float ymin = -999., float ymax = -999.) {
      if (!fRatio)
        return;
      TH1 *ratio1 = (TH1*)obj2->Clone(), *ratio2 = (TH1*)obj3->Clone();
      //ratio1->Sumw2(); ratio2->Sumw2();
      ratio1->Divide(obj1);
      ratio2->Divide(obj1);
      TCanvas::cd(2);
      ratio1->Draw("p");
      ratio2->Draw("p same");
      obj1->GetXaxis()->SetTitle("");
      if (ymin != ymax) {
        ratio1->GetYaxis()->SetRangeUser(ymin, ymax);
      }
      Prettify(ratio1);
      TCanvas::cd();
    }

    inline TH1* RatioPlot(
        TH1* obj1, const TH1* obj2 = 0, float ymin = -999., float ymax = -999., const char* plot_type = "p") {
      if (!fRatio)
        return obj1;
      TH1* ratio;
      if (obj2) {
        ratio = dynamic_cast<TH1*>(obj2->Clone());
        ratio->Divide(obj1);
      } else
        ratio = dynamic_cast<TH1*>(obj1->Clone());

      TCanvas::cd(2);
      ratio->Draw(plot_type);
      obj1->GetXaxis()->SetTitle("");
      if (ymin != ymax)
        ratio->GetYaxis()->SetRangeUser(ymin, ymax);
      Prettify(ratio);
      ratio->GetYaxis()->SetTitle("Ratio");
      TCanvas::cd();
      return ratio;
    }

    inline TGraphErrors* RatioPlot(TGraphErrors* obj1,
                                   const TGraphErrors* obj2,
                                   float ymin = -999.,
                                   float ymax = -999.) {
      if (!fRatio)
        return 0;
      TGraphErrors* ratio = new TGraphErrors;
      ratio->SetTitle(obj1->GetTitle());

      unsigned int n = 0;
      float min_x = 9.e10, max_x = -9.e10;
      for (int i = 0; i < obj1->GetN(); i++) {
        const float x1 = obj1->GetX()[i];

        for (int j = 0; j < obj2->GetN(); j++) {
          const float x2 = obj2->GetX()[j];
          if (x2 > max_x)
            max_x = x2;
          if (x2 < min_x)
            min_x = x2;

          if (fabs(x2 - x1) > 1.e-3)
            continue;
          const float y1 = obj1->GetY()[i], y1_err = obj1->GetEY()[i], y2 = obj2->GetY()[j], y2_err = obj2->GetEY()[j];
          const float y = (y2 - y1) / y1, err_y = sqrt(pow(y1_err / y1, 2) + pow(y2_err / y2, 2) * y2 / y1);
          ratio->SetPoint(n, x1, y);
          ratio->SetPointError(n, 0., err_y);
          n++;
        }
      }

      TCanvas::cd(2);
      ratio->Draw("ap");
      ratio->GetXaxis()->SetRangeUser(obj1->GetXaxis()->GetXmin(), obj1->GetXaxis()->GetXmax());
      ratio->SetMarkerStyle(20);
      if (ymin != ymax) {
        ratio->GetYaxis()->SetRangeUser(ymin, ymax);
      }
      ratio->GetXaxis()->SetLimits(min_x, max_x);
      Prettify(ratio->GetHistogram());
      obj1->GetXaxis()->SetTitle("");
      TLine l(min_x, 0., max_x, 0.);
      l.Draw();
      ratio->GetYaxis()->SetLabelSize(14);
      TCanvas::cd();

      return ratio;
    }

    inline void SetTopLabel(const char* lab = "") {
      TCanvas::cd();
      if (strcmp(lab, "") != 0)
        fTitle = lab;
      if (!fTopLabel)
        BuildTopLabel();
      else
        fTopLabel->Clear();
      fTopLabel->AddText(fTitle);
      //fTopLabel->Draw();
    }

    inline void SetLegendX1(double x) { fLegX1 = x; }
    inline void SetLegendY1(double y) { fLegY1 = y; }
    inline void AddLegendEntry(const TObject* obj, const char* title, Option_t* option = "lpf") {
      if (!fLeg)
        BuildLeg();
      fLeg->AddEntry(obj, title, option);
      const unsigned int num_entries = fLeg->GetNRows();
      if (num_entries > 3)
        fLeg->SetY1(fLeg->GetY1() - (num_entries - 3) * 0.01);
      if (num_entries > 6) {
        fLeg->SetNColumns(1 + num_entries / 6);
        fLeg->SetTextSize(0.035);
      }
    }

    inline void Save(const char* ext, const char* out_dir = ".") {
      if (strstr(ext, "pdf") == NULL)
        if (strstr(ext, "png") == NULL)
          if (strstr(ext, "root") == NULL)
            if (strstr(ext, "eps") == NULL)
              return;
      TCanvas::cd();
      if (fLeg)
        fLeg->Draw();
      if (fTopLabel)
        fTopLabel->Draw();
      TCanvas::SaveAs(Form("%s/%s.%s", out_dir, TCanvas::GetName(), ext));
    }
    inline TLegend* GetLegend() { return fLeg; }

  private:
    inline void Build() {
      TCanvas::SetLeftMargin(0.14);
      TCanvas::SetTopMargin(0.06);
      TCanvas::SetRightMargin(0.1);
      TCanvas::SetBottomMargin(0.12);
      TCanvas::SetTicks(1, 1);
      TCanvas::SetFillStyle(0);
      Pad()->SetFillStyle(0);

      SetTopLabel();
      if (fRatio)
        DivideCanvas();
    }

    inline void DivideCanvas() {
      TCanvas::Divide(1, 2);
      TPad *p1 = (TPad*)TCanvas::GetPad(1), *p2 = (TPad*)TCanvas::GetPad(2);
      p1->SetPad(0., 0.3, 1., 1.);
      p2->SetPad(0., 0.0, 1., 0.3);
      p1->SetFillStyle(0);
      p2->SetFillStyle(0);
      p1->SetLeftMargin(TCanvas::GetLeftMargin());
      p1->SetRightMargin(TCanvas::GetRightMargin());
      p2->SetLeftMargin(TCanvas::GetLeftMargin());
      p2->SetRightMargin(TCanvas::GetRightMargin());
      p1->SetTopMargin(TCanvas::GetTopMargin() + 0.025);
      p1->SetBottomMargin(0.02);
      p2->SetTopMargin(0.02);
      p2->SetBottomMargin(TCanvas::GetBottomMargin() + 0.25);
      p1->SetTicks(1, 1);
      p2->SetTicks(1, 1);
      p2->SetGrid(0, 1);
      TCanvas::cd(1);
    }

    inline void BuildTopLabel() {
      TCanvas::cd();
      fTopLabel = new TPaveText(0.5, 0.95, 0.915, 0.96, "NB NDC");
      fTopLabel->SetFillStyle(0);
      fTopLabel->SetFillColor(0);
      fTopLabel->SetLineColor(0);
      fTopLabel->SetLineStyle(0);
      fTopLabel->SetTextFont(font_type(2));
      fTopLabel->SetTextSize(0.04);
      fTopLabel->SetTextAlign(kHAlignRight + kVAlignBottom);
    }

    inline void BuildLeg() {
      if (fLeg)
        return;
      if (fRatio)
        TCanvas::cd(1);
      fLeg = new TLegend(fLegX1, fLegY1, fLegX1 + 0.45, fLegY1 + 0.15);
      fLeg->SetLineColor(kWhite);
      fLeg->SetLineWidth(0);
      fLeg->SetFillStyle(0);
      fLeg->SetTextFont(font_type(2));
      fLeg->SetTextSize(0.04);
    }
    inline float GetBinning(const TH1* h) {
      return (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin()) / h->GetXaxis()->GetNbins();
    }

    TString fTitle;
    TPaveText* fTopLabel;
    TLegend* fLeg;
    double fLegX1, fLegY1;
    bool fRatio;
  };
}  // namespace cepgen

#endif
