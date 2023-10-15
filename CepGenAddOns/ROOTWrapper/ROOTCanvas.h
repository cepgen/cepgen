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
#include <TMath.h>
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
    inline explicit ROOTPaveText(float x1, float y1, float x2, float y2, const std::string& text = "")
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
        : TCanvas(name.c_str(), "", 600, 600), ratio_(ratio) {
      gStyle->SetOptStat(0);
      gStyle->SetGridColor(17);
      gStyle->SetEndErrorSize(0);
      SetTopLabel(title);
      Build();
    }

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

    inline void Prettify(THStack* hs) {
      Prettify(hs->GetHistogram());
      if (!ratio_)
        return;
      auto* objarr = hs->GetHists();
      if (objarr->GetEntries() < 2)
        return;
      TH1* denom = nullptr;
      std::vector<TH1*> numers{};
      for (int i = 0; i < objarr->GetEntries(); ++i)
        if (i == 0) {  // reference is conventionally the first histogram
          denom = dynamic_cast<TH1*>(objarr->At(i)->Clone());
          denom->GetXaxis()->SetTitle(hs->GetHistogram()->GetXaxis()->GetTitle());
        } else
          numers.emplace_back(dynamic_cast<TH1*>(objarr->At(i)->Clone()));
      RatioPlot(denom, numers);
    }
    inline void Prettify(TMultiGraph* mg) {
      Prettify(mg->GetHistogram());
      if (!ratio_)
        return;
      auto* list = mg->GetListOfGraphs();
      if (list->GetEntries() < 2)
        return;
      TGraphErrors* denom = nullptr;
      std::vector<TGraphErrors*> numers{};
      double x_min{1.e10}, x_max{-1.e10};
      for (int i = 0; i < list->GetEntries(); ++i) {
        TGraphErrors* gre{nullptr};
        if (strcmp(list->At(i)->ClassName(), "TGraph") == 0) {
          auto* gr = dynamic_cast<TGraph*>(list->At(i));
          gre = new TGraphErrors(gr->GetN(), gr->GetX(), gr->GetY());
          gre->SetLineColor(gr->GetLineColor());
          gre->SetLineWidth(gr->GetLineWidth());
          gre->SetLineStyle(gr->GetLineStyle());
          gre->SetTitle(gr->GetTitle());
        } else if (strcmp(list->At(i)->ClassName(), "TGraphErrors") == 0)
          gre = dynamic_cast<TGraphErrors*>(list->At(i));
        gre->SetTitle(mg->GetTitle());
        x_min = TMath::Min(TMath::MinElement(gre->GetN(), gre->GetX()), x_min);
        x_max = TMath::Max(TMath::MaxElement(gre->GetN(), gre->GetX()), x_max);
        if (i == 0) {  // reference is conventionally the first graph
          denom = gre;
        } else
          numers.emplace_back(gre);
      }
      RatioPlot(denom, numers, x_min, x_max);
      mg->GetXaxis()->SetRangeUser(x_min, x_max);
    }

    inline std::vector<TH1*> RatioPlot(TH1* denom,
                                       const std::vector<TH1*>& numers,
                                       float xmin = -999.,
                                       float xmax = -999.,
                                       float ymin = -999.,
                                       float ymax = -999.,
                                       Option_t* draw_style = "hist") {
      std::vector<TH1*> ratios{};
      if (!ratio_)
        return ratios;
      TCanvas::cd(2);
      auto* hs = Make<THStack>();  // garbage collected
      for (const auto& numer : numers) {
        auto* ratio = dynamic_cast<TH1*>(numer->Clone("ratio"));
        ratio->Divide(denom);
        auto* ratio_shadow = dynamic_cast<TH1*>(ratio->Clone("ratio_shadow"));
        ratio_shadow->SetFillColorAlpha(ratio->GetLineColor(), 0.25);
        hs->Add(ratio_shadow, "e2");
        hs->Add(ratio, draw_style);
        ratios.emplace_back(ratio);
      }
      pads_.at(1)->SetLogy(false);
      hs->Draw("nostack");
      if (xmin == xmax) {
        xmin = denom->GetXaxis()->GetXmin();
        xmax = denom->GetXaxis()->GetXmax();
      }
      TLine l;
      l.SetLineWidth(2);
      l.DrawLine(xmin, 1., xmax, 1.);
      auto* hst = hs->GetHistogram();
      Prettify(hst);
      hst->GetXaxis()->SetTitle(denom->GetXaxis()->GetTitle());
      hst->GetXaxis()->SetTitleOffset(0.);
      hst->GetXaxis()->SetTickSize(0.065);
      hst->GetXaxis()->SetRangeUser(xmin, xmax);
      hst->GetYaxis()->SetTitle("Ratio");
      if (ymin != ymax)
        hst->GetYaxis()->SetRangeUser(ymin, ymax);
      else
        hst->GetYaxis()->SetRangeUser(TMath::Max(-0.1, hst->GetYaxis()->GetXmin()),
                                      TMath::Min(2.1, hst->GetYaxis()->GetXmax()));
      denom->GetXaxis()->SetTitle("");
      TCanvas::cd(1);
      return ratios;
    }

    inline std::vector<TGraphErrors*> RatioPlot(TGraphErrors* denom,
                                                const std::vector<TGraphErrors*>& numers,
                                                float xmin = -999.,
                                                float xmax = -999.,
                                                float ymin = -999.,
                                                float ymax = -999.) {
      std::vector<TGraphErrors*> ratios{};
      if (!ratio_)
        return ratios;
      auto* mg = Make<TMultiGraph>();
      const auto *xd = denom->GetX(), *yd = denom->GetY(), *yde = denom->GetEY();
      for (const auto& numer : numers) {
        if (numer->GetN() != denom->GetN())
          continue;
        const auto *xn = numer->GetX(), *yn = numer->GetY(), *yne = numer->GetEY();
        auto* ratio = new TGraphErrors();
        ratio->SetTitle(denom->GetTitle());
        for (int i = 0; i < denom->GetN(); i++) {
          const float xd_val = xd[i], yd_val = yd[i], yd_err = yde[i];
          for (int j = 0; j < numer->GetN(); ++j) {
            const float xn_val = xn[j], yn_val = yn[j], yn_err = yne[j];
            if ((xn_val == 0. && xd_val == 0.) || fabs(1. - xd_val / xn_val) < 1.e-2) {
              if (yd_val == 0. || yn_val == 0.)
                break;
              const float y = yn_val / yd_val, err_y = std::hypot(yn_err / yn_val, yd_err / yd_val) * y;
              const auto n = ratio->GetN();
              ratio->SetPoint(n, xd_val, y);
              ratio->SetPointError(n, 0., err_y);
              break;
            }
          }
        }
        mg->Add(ratio);
        ratio->SetLineColor(numer->GetLineColor());
        ratio->SetLineWidth(numer->GetLineWidth());
        ratio->SetLineStyle(numer->GetLineStyle());
        ratios.emplace_back(ratio);
      }
      TCanvas::cd(2);
      mg->Draw("al");
      Prettify(mg->GetHistogram());
      if (xmin != xmax)
        mg->GetXaxis()->SetRangeUser(xmin, xmax);
      else
        mg->GetXaxis()->SetRangeUser(denom->GetXaxis()->GetXmin(), denom->GetXaxis()->GetXmax());
      mg->GetXaxis()->SetTitle(denom->GetXaxis()->GetTitle());
      mg->GetXaxis()->SetTitleOffset(0.);
      mg->GetXaxis()->SetTickSize(0.065);
      if (ymin != ymax)
        mg->GetYaxis()->SetRangeUser(ymin, ymax);
      else
        mg->GetYaxis()->SetRangeUser(TMath::Max(-0.1, mg->GetYaxis()->GetXmin()),
                                     TMath::Min(2.1, mg->GetYaxis()->GetXmax()));
      mg->GetYaxis()->SetTitle("Ratio");
      denom->GetXaxis()->SetTitle("");
      TLine l(denom->GetXaxis()->GetXmin(), 1., denom->GetXaxis()->GetXmax(), 1.);
      l.Draw();
      TCanvas::cd(1);
      return ratios;
    }
    /// Specify the text to show on top of the canvas
    inline void SetTopLabel(const std::string& lab) {
      TCanvas::cd();
      std::string title = "CepGen v" + version::tag;
      if (!lab.empty())
        title += " - " + lab;
      if (!top_label_)
        BuildTopLabel();
      else
        top_label_->Clear();
      top_label_->AddText(title.data());
      //top_label_->Draw();
    }
    inline void SetGrid(int x = true, int y = true) { pads_.at(0)->SetGrid(x, y); }
    inline void SetLogx(int log = true) {
      for (auto& pad : pads_)
        pad->SetLogx(log);
    }
    inline void SetLogy(int log = true) { pads_.at(0)->SetLogy(log); }
    inline void SetLogz(int log = true) { pads_.at(0)->SetLogz(log); }

    /// Set the placement strategy for the legend
    inline void SetLegendMode(const std::string& mode) { leg_mode_ = mode; }
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
        leg_width_ = 0.55;
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
        if (TPad::PlaceBox(leg_.get(), leg_width_ * 1.15, leg_height_, leg_x1_, leg_y1_, leg_mode_.data())) {
          leg_y1_ = std::min(leg_y1_, 0.9 - leg_height_);
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
    inline void Place(TLegend* leg, const Option_t* mode = "lt") {
      if (!leg)
        return;
      double leg_x, leg_y;
      const auto leg_width = leg->GetX2() - leg->GetX1(), leg_height = leg->GetY2() - leg->GetY1();
      if (TPad::PlaceBox(leg, leg_width * 1.15, leg_height, leg_x, leg_y, mode)) {
        leg->SetX1(leg_x);
        leg->SetX2(leg_x + leg_width);
        leg->SetY1(leg_y);
        leg->SetY2(leg_y + leg_height);
      }
      leg->Draw();
    }
    /// Garbage collector-like TObjects producer
    template <typename T, typename... Args>
    inline T* Make(Args&&... args) {
      grb_obj_.emplace_back(new T(std::forward<Args>(args)...));
      return dynamic_cast<T*>(grb_obj_.rbegin()->get());
    }

  private:
    /// Prepare the canvas for later drawing
    inline void Build() {
      TCanvas::SetLeftMargin(0.14);
      TCanvas::SetTopMargin(0.06);
      TCanvas::SetRightMargin(0.1);
      TCanvas::SetBottomMargin(0.12);
      TCanvas::SetTicks(1, 1);
      TCanvas::SetFillStyle(0);
      TCanvas::Pad()->SetFillStyle(0);
      if (ratio_)
        DivideCanvas();
    }
    /// Divide the canvas into two sub-pads if a ratio plot is to be shown
    inline void DivideCanvas() {
      TCanvas::Pad()->Divide(1, 2);
      pads_.clear();
      // main pad
      auto* p1 = dynamic_cast<TPad*>(TCanvas::GetPad(1));
      p1->SetPad(0., 0.3, 1., 1.);
      p1->SetFillStyle(0);
      p1->SetLeftMargin(TCanvas::GetLeftMargin());
      p1->SetRightMargin(TCanvas::GetRightMargin());
      p1->SetTopMargin(TCanvas::GetTopMargin() + 0.025);
      p1->SetBottomMargin(0.02);
      p1->SetTicks(1, 1);
      pads_.emplace_back(p1);
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
      pads_.emplace_back(p2);
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

    const bool ratio_;
    std::string leg_mode_{"rt"};
    double leg_x1_{0.5}, leg_y1_{0.75};
    double leg_width_{0.45}, leg_height_{0.15};
    std::unique_ptr<TLegend> leg_;
    std::unique_ptr<ROOTPaveText> top_label_;
    std::vector<std::unique_ptr<TObject> > grb_obj_;
    std::vector<TPad*> pads_{};
  };
  const std::vector<int> ROOTCanvas::colours = {
      kBlack, kRed + 1, kBlue - 2, kGreen + 1, kOrange + 1, kAzure + 1, kMagenta + 1, kCyan + 3, kPink + 5};
}  // namespace cepgen

#endif
