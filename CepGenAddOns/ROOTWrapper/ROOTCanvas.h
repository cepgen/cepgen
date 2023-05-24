/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#ifndef CepGenAddOns_ROOTWrapper_ROOTCanvas_h
#define CepGenAddOns_ROOTWrapper_ROOTCanvas_h

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <cstring>
#include <vector>

#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  /// A "prettified" text box object
  class ROOTPaveText : public TPaveText {
  public:
    inline ROOTPaveText(float x1, float y1, float x2, float y2, const std::string& text = "")
        : TPaveText(x1, y1, x2, y2, "NB NDC") {
      TPaveText::SetTextAlign(kHAlignLeft + kVAlignTop);
      if (!text.empty()) {
        TString txt = text;
        if (txt.Contains("\\")) {
          TObjArray* tok = txt.Tokenize("\\");
          for (int i = 0; i < tok->GetEntries(); ++i)
            TPaveText::AddText(dynamic_cast<TObjString*>(tok->At(i))->String());
        } else
          TPaveText::AddText(txt);
      }
      TPaveText::SetFillColor(0);
      TPaveText::SetFillStyle(0);
      TPaveText::SetLineColor(0);
      TPaveText::SetLineWidth(0);
      TPaveText::SetShadowColor(0);
      TPaveText::SetTextFont(fontType(2));
      TPaveText::SetTextSize(0.058);
    }

    /// Force font to be Times New Roman-style
    inline static int fontType(int mode) { return 130 + mode; }
  };

  /// A "prettified" generic figure canvas
  class ROOTCanvas : public TCanvas {
  public:
    static const std::vector<int> colours;

    /// Build a canvas from its name, title, and attributes
    /// \param[in] name Canvas name (and subsequently filename on save)
    /// \param[in] ratio Divide the canvas into a main and ratio plots subparts?
    explicit inline ROOTCanvas(const std::string& name, const std::string& title = "", bool ratio = false)
        : TCanvas(name.c_str(), "", 600, 600), title_(title), ratio_(ratio) {
      gStyle->SetOptStat(0);
      gStyle->SetGridColor(17);
      gStyle->SetEndErrorSize(0);
      Build();
    }
    inline ~ROOTCanvas() {}

    /// Set horizontal canvas width
    inline void SetSize(float size = 600) { TCanvas::SetCanvasSize(size, 600); }

    /// Draw main plot attributes in a pretty manner
    inline void Prettify(TH1* obj) {
      auto* x = dynamic_cast<TAxis*>(obj->GetXaxis());
      x->CenterTitle();
      x->SetLabelFont(ROOTPaveText::fontType(3));
      x->SetLabelSize(20);
      x->SetTitleFont(ROOTPaveText::fontType(3));
      x->SetTitleSize(29);
      if (ratio_) {
        x->SetTitleOffset(2.5);
        x->SetLabelOffset(0.02);
      }
      x->SetTickLength(0.03);
      auto* y = dynamic_cast<TAxis*>(obj->GetYaxis());
      y->CenterTitle();
      y->SetLabelFont(ROOTPaveText::fontType(3));
      y->SetLabelSize(20);
      y->SetTitleFont(ROOTPaveText::fontType(3));
      y->SetTitleSize(29);
      y->SetTitleOffset(1.3);
      y->SetTickLength(0.03);
      auto* z = dynamic_cast<TAxis*>(obj->GetZaxis());
      z->CenterTitle();
      z->SetLabelFont(ROOTPaveText::fontType(3));
      z->SetLabelSize(16);
      z->SetTitleFont(ROOTPaveText::fontType(3));
      z->SetTitleSize(29);

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
        x->SetTitle(x_title);
        y->SetTitle(y_title);
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

    inline std::vector<TH1*> RatioPlot(TH1* denom,
                                       const std::vector<TH1*>& numers,
                                       float ymin = -999.,
                                       float ymax = -999.,
                                       Option_t* draw_style = "hist") {
      std::vector<TH1*> ratios;
      if (!ratio_)
        return ratios;
      auto* hs = Make<THStack>();  // garbage collected
      for (const auto& numer : numers) {
        auto* ratio = dynamic_cast<TH1*>(numer->Clone("ratio"));
        //ratio->Divide(denom);
        //ratio->Sumw2();
        hs->Add(ratio, draw_style);
        ratios.emplace_back(ratio);
      }
      TCanvas::cd(2);
      hs->Draw("nostack");
      auto* hst = hs->GetHistogram();
      Prettify(hst);
      if (ymin != ymax)
        hst->GetYaxis()->SetRangeUser(ymin, ymax);
      hst->GetYaxis()->SetTitle("Ratio");
      printf("%s::%s\n", denom->GetTitle(), denom->GetXaxis()->GetTitle());
      hst->GetXaxis()->SetTitle(denom->GetXaxis()->GetTitle());
      const double xmin = denom->GetXaxis()->GetXmin(), xmax = denom->GetXaxis()->GetXmax();
      hst->GetXaxis()->SetLimits(xmin, xmax);
      TLine l;
      l.SetLineWidth(2);
      l.DrawLine(xmin, 1., xmax, 1.);
      denom->GetXaxis()->SetTitle("");
      TCanvas::cd();
      return ratios;
    }

    inline TGraphErrors* RatioPlot(TGraphErrors* obj1,
                                   const TGraphErrors* obj2,
                                   float ymin = -999.,
                                   float ymax = -999.) {
      if (!ratio_)
        return 0;
      auto* ratio = Make<TGraphErrors>();
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
    /// Specify the text to show on top of the canvas
    inline void SetTopLabel(const std::string& lab) {
      TCanvas::cd();
      title_ = "CepGen v" + version::tag;
      if (!lab.empty())
        title_ += " - " + lab;
      if (!top_label_)
        BuildTopLabel();
      else
        top_label_->Clear();
      top_label_->AddText(title_.data());
      //top_label_->Draw();
    }

    /// Set the horizontal coordinate of the low-left part of the legend object
    /// \note To be called before the first legend entry is added
    inline void SetLegendX1(double x) {
      if (leg_)
        perror("SetLegendX1");
      leg_x1_ = x;
    }
    /// Set the vertical coordinate of the low-left part of the legend object
    /// \note To be called before the first legend entry is added
    inline void SetLegendY1(double y) {
      if (leg_)
        perror("SetLegendY1");
      leg_y1_ = y;
    }
    /// Add one new entry to the legend object
    inline void AddLegendEntry(const TObject* obj, const std::string& title, Option_t* option = "lpf") {
      if (!leg_)
        BuildLeg();
      leg_->AddEntry(obj, title.c_str(), option);
      const unsigned int num_entries = leg_->GetNRows();
      if (num_entries > 3)
        leg_->SetY1(leg_->GetY1() - (num_entries - 3) * 0.01);
      if (num_entries > 6) {
        leg_->SetNColumns(1 + num_entries / 6);
        leg_->SetTextSize(0.035);
      }
    }
    /// Save the canvas in an external file
    inline void Save(const std::string& ext, const std::string& out_dir = ".") {
      const auto& extensions = utils::split(ext, ',');
      if (extensions.empty())
        return;
      TCanvas::cd();
      if (top_label_)
        top_label_->Draw();
      if (leg_) {
        double leg_x1, leg_y1;
        if (TPad::PlaceBox(leg_.get(), leg_width_, leg_height_, leg_x1, leg_y1, "lb")) {
          leg_x1_ = std::max(leg_x1_, leg_x1);
          leg_y1_ = std::max(std::min(leg_y1_, leg_y1), 0.9 - leg_height_);
          leg_->SetX1(leg_x1_);
          leg_->SetX2(leg_x1_ + leg_width_);
          leg_->SetY1(leg_y1_);
          leg_->SetY2(leg_y1_ + leg_height_);
        }
        leg_->Draw();
      }
      for (const auto& extension : extensions)
        TCanvas::SaveAs(Form("%s/%s.%s", out_dir.c_str(), TCanvas::GetName(), extension.c_str()));
    }
    /// Retrieve the legend object (if produced)
    inline TLegend* GetLegend() { return leg_.get(); }

  private:
    /// Prepare the canvas for later drawing
    inline void Build() {
      TCanvas::SetLeftMargin(0.14);
      TCanvas::SetTopMargin(0.06);
      TCanvas::SetRightMargin(0.1);
      TCanvas::SetBottomMargin(0.12);
      TCanvas::SetTicks(1, 1);
      TCanvas::SetFillStyle(0);
      Pad()->SetFillStyle(0);

      SetTopLabel("");
      if (ratio_)
        DivideCanvas();
    }
    /// Divide the canvas into two sub-pads if a ratio plot is to be shown
    inline void DivideCanvas() {
      TCanvas::Divide(1, 2);
      // main pad
      auto* p1 = dynamic_cast<TPad*>(TCanvas::GetPad(1));
      p1->SetPad(0., 0.3, 1., 1.);
      p1->SetFillStyle(0);
      p1->SetLeftMargin(TCanvas::GetLeftMargin());
      p1->SetRightMargin(TCanvas::GetRightMargin());
      p1->SetTopMargin(TCanvas::GetTopMargin() + 0.025);
      p1->SetBottomMargin(0.02);
      p1->SetTicks(1, 1);
      // ratio plot(s) pad
      auto* p2 = dynamic_cast<TPad*>(TCanvas::GetPad(2));
      p2->SetPad(0., 0.0, 1., 0.3);
      p2->SetFillStyle(0);
      p2->SetLeftMargin(TCanvas::GetLeftMargin());
      p2->SetRightMargin(TCanvas::GetRightMargin());
      p2->SetTopMargin(0.02);
      p2->SetBottomMargin(TCanvas::GetBottomMargin() + 0.25);
      p2->SetTicks(1, 1);
      p2->SetGrid(0, 1);
      // roll back to main pad
      TCanvas::cd(1);
    }
    /// Build the text box on top of the canvas
    inline void BuildTopLabel() {
      TCanvas::cd();
      top_label_.reset(new ROOTPaveText(0.5, 0.95, 0.915, 0.96));
      top_label_->SetTextSize(0.04);
      top_label_->SetTextAlign(kHAlignRight + kVAlignBottom);
    }
    /// Build the legend object if not already done
    inline void BuildLeg() {
      if (leg_)
        return;
      if (ratio_)
        TCanvas::cd(1);
      leg_.reset(new TLegend);
      leg_->SetLineColor(kWhite);
      leg_->SetLineWidth(0);
      leg_->SetFillStyle(0);
      leg_->SetTextFont(ROOTPaveText::fontType(2));
      leg_->SetTextSize(0.04);
    }
    /// Retrieve the bin size for a histogram
    inline float GetBinning(const TH1* hist) {
      return (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / hist->GetXaxis()->GetNbins();
    }
    /// Garbage collector-like TObjects producer
    template <typename T, typename... Args>
    inline T* Make(Args&&... args) {
      grb_obj_.emplace_back(new T(std::forward<Args>(args)...));
      return dynamic_cast<T*>(grb_obj_.rbegin()->get());
    }

    std::string title_;
    bool ratio_;
    double leg_x1_{0.5}, leg_y1_{0.75};
    double leg_width_{0.45}, leg_height_{0.15};
    std::unique_ptr<TLegend> leg_;
    std::unique_ptr<ROOTPaveText> top_label_;
    std::vector<std::unique_ptr<TObject> > grb_obj_;
  };
  const std::vector<int> ROOTCanvas::colours = {
      kBlack, kRed + 1, kBlue - 2, kGreen + 1, kOrange + 1, kAzure + 1, kMagenta + 1, kCyan + 3, kPink + 5};
}  // namespace cepgen

#endif
