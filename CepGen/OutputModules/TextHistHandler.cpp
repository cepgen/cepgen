/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class TextHistHandler : public ExportModule {
    public:
      explicit TextHistHandler(const ParametersList&);
      ~TextHistHandler();

      static ParametersDescription description();

      void initialise(const Parameters&) override;
      void setCrossSection(double cross_section, double) override { cross_section_ = cross_section; }
      void operator<<(const Event&) override;

    private:
      std::ofstream file_;
      const std::unique_ptr<utils::Drawer> drawer_;
      //--- variables definition
      const bool show_hists_, save_hists_;
      const std::string filename_;

      const utils::EventBrowser browser_;

      double cross_section_{1.};
      unsigned long num_evts_{0ul};

      /// centre-of-mass energy
      double sqrts_{0.};

      /// 1D histogram definition
      struct Hist1DInfo {
        std::string var;
        utils::Hist1D hist;
        bool log;
      };
      /// List of 1D histograms
      std::vector<Hist1DInfo> hists_;
      /// 2D histogram definition
      struct Hist2DInfo {
        std::string var1, var2;
        utils::Hist2D hist;
        bool log;
      };
      /// List of 2D histograms
      std::vector<Hist2DInfo> hists2d_;
    };

    TextHistHandler::TextHistHandler(const ParametersList& params)
        : ExportModule(params),
          drawer_(utils::DrawerFactory::get().build("text", params)),
          show_hists_(steer<bool>("showHistograms")),
          save_hists_(steer<bool>("saveHistograms")),
          filename_(steer<std::string>("filename")) {
      //--- extract list of variables to be plotted in histogram
      const auto& hist_vars = steer<ParametersList>("histVariables");
      for (const auto& key : hist_vars.keys()) {
        const auto& vars = utils::split(key, ':');
        if (vars.size() < 1 || vars.size() > 2)
          throw CG_FATAL("TextHistHandler") << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvar = hist_vars.get<ParametersList>(key);
        const auto& log = hvar.get<bool>("log");
        auto name = utils::replace_all(utils::replace_all(key, ")", ""), "(", "_");
        if (vars.size() == 1) {  // 1D histogram
          const auto& xbins = hvar.get<std::vector<double> >("xbins");
          const auto title = "d(sig)/d" + vars.at(0) + " (pb/bin)";
          if (xbins.size() > 1)
            hists_.emplace_back(Hist1DInfo{vars.at(0), utils::Hist1D(xbins, name, title), log});
          else if (hvar.get<Limits>("xrange").valid()) {
            const auto& nbins = (hvar.get<int>("nbins") > 0 ? hvar.get<int>("nbins") : hvar.get<int>("nbinsX"));
            hists_.emplace_back(
                Hist1DInfo{vars.at(0), utils::Hist1D(nbins, hvar.get<Limits>("xrange"), name, title), log});
          } else {
            CG_WARNING("TextHistHandler")
                << "Neither xrange nor xbins found in parameters for 1D plot of variable \"" << vars.at(0) << "\".";
            continue;
          }
          auto& hist = hists_.rbegin()->hist;
          hist.setXlabel(vars.at(0));
          hist.setYlabel("d(sig)/d" + vars.at(0) + " (pb/bin)");
        } else if (vars.size() == 2) {  // 2D histogram
          const auto &xbins = hvar.get<std::vector<double> >("xbins"), &ybins = hvar.get<std::vector<double> >("ybins");
          name = utils::replace_all(name, ":", "_");
          const auto title = "d^2(sig)/d" + vars.at(0) + "/d" + vars.at(1) + " (pb/bin)";
          if (xbins.size() > 1 && ybins.size() > 1)
            hists2d_.emplace_back(Hist2DInfo{vars.at(0), vars.at(1), utils::Hist2D(xbins, ybins, name, title), log});
          else if (hvar.get<Limits>("xrange").valid()) {
            const auto& nbinsx = (hvar.get<int>("nbins") > 0 ? hvar.get<int>("nbins") : hvar.get<int>("nbinsX"));
            hists2d_.emplace_back(Hist2DInfo{vars.at(0),
                                             vars.at(1),
                                             utils::Hist2D(nbinsx,
                                                           hvar.get<Limits>("xrange"),
                                                           hvar.get<int>("nbinsY"),
                                                           hvar.get<Limits>("yrange"),
                                                           name,
                                                           title),
                                             log});
          } else {
            CG_WARNING("TextHistHandler")
                << "Neither (x/y)range nor (x/y)bins found in parameters for 1D plot of variables \"" << vars << "\".";
            continue;
          }
          auto& hist = hists2d_.rbegin()->hist;
          hist.setXlabel(vars.at(0));
          hist.setYlabel(vars.at(1));
        }
      }
      if (save_hists_ && !hists_.empty())
        file_.open(filename_);
    }

    TextHistHandler::~TextHistHandler() {
      //--- histograms printout
      if (!show_hists_ && !save_hists_)
        return;
      for (auto& h_var : hists_) {
        h_var.hist.scale(cross_section_ / (num_evts_ + 1));
        std::ostringstream os;
        drawer_->draw(h_var.hist, h_var.log ? utils::Drawer::Mode::logy : utils::Drawer::Mode::none);
        if (show_hists_)
          CG_INFO("TextHistHandler") << os.str();
        if (save_hists_)
          file_ << "\n" << os.str() << "\n";
      }
      for (const auto& h_var : hists2d_) {
        std::ostringstream os;
        drawer_->draw(h_var.hist, h_var.log ? utils::Drawer::Mode::logy : utils::Drawer::Mode::none);
        if (show_hists_)
          CG_INFO("TextHistHandler") << os.str();
        if (save_hists_)
          file_ << "\n" << os.str() << "\n";
      }
      if (save_hists_)
        CG_INFO("TextHistHandler") << "Saved " << utils::s("histogram", hists_.size(), true) << " into \"" << filename_
                                   << "\".";
    }

    void TextHistHandler::initialise(const Parameters& params) {
      sqrts_ = params.kinematics().incomingBeams().sqrtS();
      num_evts_ = 0ul;
      if (save_hists_ && !hists_.empty())
        file_ << banner(params, "#") << "\n";
    }

    void TextHistHandler::operator<<(const Event& ev) {
      //--- increment the corresponding histograms
      for (auto& h_var : hists_)
        h_var.hist.fill(browser_.get(ev, h_var.var));
      for (auto& h_var : hists2d_)
        h_var.hist.fill(browser_.get(ev, h_var.var1), browser_.get(ev, h_var.var2));
      ++num_evts_;
    }

    ParametersDescription TextHistHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("Text-based histogramming tool");
      desc.add<std::string>("filename", "output.hists.txt").setDescription("Output filename for histogram dump");
      desc.add<bool>("showHistograms", true).setDescription("Show the histogram(s) at the end of the run?");
      desc.add<bool>("saveHistograms", false).setDescription("Save the histogram(s) at the end of the run?");
      // per-histogram default parameters
      ParametersDescription hist_desc;
      // x-axis attributes
      hist_desc.add<std::vector<double> >("xbins", {}).setDescription("x-axis bins definition");
      hist_desc.add<int>("nbins", 25).setDescription("Bins multiplicity for x-axis");
      hist_desc.add<int>("nbinsX", -1).setDescription("Bins multiplicity for x-axis");
      hist_desc.add<Limits>("xrange", Limits{0., 1.}).setDescription("Minimum-maximum range for x-axis");
      // y-axis attributes
      hist_desc.add<std::vector<double> >("ybins", {}).setDescription("y-axis bins definition");
      hist_desc.add<int>("nbinsY", 50).setDescription("Bins multiplicity for y-axis");
      hist_desc.add<Limits>("yrange", Limits{0., 1.}).setDescription("Minimum-maximum range for y-axis");
      hist_desc.add<bool>("log", false).setDescription("Plot logarithmic axis?");
      desc.addParametersDescriptionVector("histVariables", hist_desc, {})
          .setDescription("Histogram definition for 1/2 variable(s)");
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("text", TextHistHandler)
