/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/Plotter.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class TextHandler : public ExportModule {
    public:
      explicit TextHandler(const ParametersList&);
      ~TextHandler();

      static ParametersDescription description();

      void initialise(const Parameters&) override;
      void setCrossSection(double cross_section, double) override { cross_section_ = cross_section; }
      void operator<<(const Event&) override;

    private:
      std::ofstream file_, hist_file_;
      std::string hist_filename_;
      //--- variables definition
      const std::vector<std::string> variables_;
      const bool save_banner_, save_variables_;
      const bool show_hists_, save_hists_;
      const std::string separator_;

      const utils::EventBrowser browser_;

      std::ostringstream oss_vars_;

      double cross_section_{1.};

      //--- kinematic variables
      double sqrts_{0.};
      unsigned long num_evts_{0ul};
      struct Hist1DInfo {
        std::string var;
        utils::Hist1D hist;
      };
      std::vector<Hist1DInfo> hists_;
      struct Hist2DInfo {
        std::string var1, var2;
        utils::Hist2D hist;
      };
      std::vector<Hist2DInfo> hists2d_;
    };

    TextHandler::TextHandler(const ParametersList& params)
        : ExportModule(params),
          file_(params.get<std::string>("filename")),
          hist_filename_(params.get<std::string>("histFilename")),
          variables_(params.get<std::vector<std::string> >("variables")),
          save_banner_(params.get<bool>("saveBanner")),
          save_variables_(params.get<bool>("saveVariables")),
          show_hists_(params.get<bool>("showHistograms")),
          save_hists_(params.get<bool>("saveHistograms")),
          separator_(params.get<std::string>("separator")) {
      //--- first extract list of variables to store in output file
      oss_vars_.clear();
      std::string sep;
      for (const auto& var : variables_)
        oss_vars_ << sep << var, sep = separator_;
      //--- then extract list of variables to be plotted in histogram
      const auto& hist_vars = params.get<ParametersList>("histVariables");
      for (const auto& key : hist_vars.keys()) {
        const auto& vars = utils::split(key, ':');
        if (vars.size() < 1 || vars.size() > 2)
          throw CG_FATAL("TextHandler") << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvar = hist_vars.get<ParametersList>(key);
        if (vars.size() == 1) {  // 1D histogram
          if (hvar.has<std::vector<double> >("xbins"))
            hists_.emplace_back(Hist1DInfo{vars.at(0), utils::Hist1D(hvar.get<std::vector<double> >("xbins"))});
          else if (hvar.has<Limits>("xrange")) {
            const auto& nbins = (hvar.get<int>("nbins") > 0 ? hvar.get<int>("nbins") : hvar.get<int>("nbinsX"));
            hists_.emplace_back(Hist1DInfo{vars.at(0), utils::Hist1D(nbins, hvar.get<Limits>("xrange"))});
          } else {
            CG_WARNING("TextHandler") << "Neither xrange nor xbins found in parameters for 1D plot of variable \""
                                      << vars.at(0) << "\".";
            continue;
          }
          auto& hist = hists_.rbegin()->hist;
          hist.setLog(hvar.get<bool>("log"));
          hist.setName(key);
          hist.setXlabel(vars.at(0));
          hist.setYlabel("d(sig)/d" + vars.at(0) + " (pb/bin)");
        } else if (vars.size() == 2) {  // 2D histogram
          if (hvar.has<std::vector<double> >("xbins") && hvar.has<std::vector<double> >("ybins"))
            hists2d_.emplace_back(Hist2DInfo{
                vars.at(0),
                vars.at(1),
                utils::Hist2D(hvar.get<std::vector<double> >("xbins"), hvar.get<std::vector<double> >("ybins"))});
          else if (hvar.has<Limits>("xrange")) {
            const auto& nbinsx = (hvar.get<int>("nbins") > 0 ? hvar.get<int>("nbins") : hvar.get<int>("nbinsX"));
            CG_WARNING("") << nbinsx << ": " << hvar;
            hists2d_.emplace_back(Hist2DInfo{
                vars.at(0),
                vars.at(1),
                utils::Hist2D(
                    nbinsx, hvar.get<Limits>("xrange"), hvar.get<int>("nbinsY"), hvar.get<Limits>("yrange"))});
          } else {
            CG_WARNING("TextHandler")
                << "Neither (x/y)range nor (x/y)bins found in parameters for 1D plot of variables \"" << vars << "\".";
            continue;
          }
          auto& hist = hists2d_.rbegin()->hist;
          hist.setName(key);
          hist.setXlabel(vars.at(0));
          hist.setYlabel(vars.at(1));
          hist.setName("d^2(sig)/d" + vars.at(0) + "/d" + vars.at(1) + " (pb/bin)");
          hist.setLog(hvar.get<bool>("log"));
        }
      }
      if (save_hists_ && !hists_.empty())
        hist_file_.open(hist_filename_);
    }

    TextHandler::~TextHandler() {
      //--- finalisation of the output file
      file_.close();
      //--- histograms printout
      if (!show_hists_ && !save_hists_)
        return;
      for (auto& h_var : hists_) {
        h_var.hist.scale(cross_section_ / (num_evts_ + 1));
        std::ostringstream os;
        h_var.hist.draw(os);
        if (show_hists_)
          CG_INFO("TextHandler") << os.str();
        if (save_hists_)
          hist_file_ << "\n" << os.str() << "\n";
      }
      for (const auto& h_var : hists2d_) {
        std::ostringstream os;
        h_var.hist.draw(os);
        if (show_hists_)
          CG_INFO("TextHandler") << os.str();
        if (save_hists_)
          hist_file_ << "\n" << os.str() << "\n";
      }
      if (save_hists_)
        CG_INFO("TextHandler") << "Saved " << utils::s("histogram", hists_.size()) << " into \"" << hist_filename_
                               << "\".";
    }

    void TextHandler::initialise(const Parameters& params) {
      sqrts_ = params.kinematics().incomingBeams().sqrtS();
      num_evts_ = 0ul;
      if (save_banner_)
        file_ << banner(params, "#") << "\n";
      if (save_variables_)
        file_ << "# " << oss_vars_.str() << "\n";
      if (save_hists_ && !hists_.empty())
        hist_file_ << banner(params, "#") << "\n";
    }

    void TextHandler::operator<<(const Event& ev) {
      //--- write down the variables list in the file
      if (!variables_.empty()) {
        std::string sep;
        for (const auto& var : variables_)
          file_ << sep << browser_.get(ev, var), sep = separator_;
        file_ << "\n";
      }
      //--- increment the corresponding histograms
      for (auto& h_var : hists_)
        h_var.hist.fill(browser_.get(ev, h_var.var));
      for (auto& h_var : hists2d_)
        h_var.hist.fill(browser_.get(ev, h_var.var1), browser_.get(ev, h_var.var2));
      ++num_evts_;
    }

    ParametersDescription TextHandler::description() {
      auto desc = ExportModule::description();
      desc.setDescription("Text-based histogramming tool");
      desc.add<std::string>("filename", "output.txt").setDescription("Output filename for variables dump");
      desc.add<std::string>("histFilename", "output.hists.txt").setDescription("Output filename for histogram dump");
      desc.add<std::vector<std::string> >("variables", {}).setDescription("List of variables to dump");
      desc.add<bool>("saveBanner", true).setDescription("Also save the boilerplate in output files?");
      desc.add<bool>("saveVariables", true).setDescription("Save the variable(s) into an output file?");
      desc.add<bool>("showHistograms", true).setDescription("Show the histogram(s) at the end of the run?");
      desc.add<bool>("saveHistograms", false).setDescription("Save the histogram(s) at the end of the run?");
      desc.add<std::string>("separator", "\t").setDescription("Base separator in output file");
      // per-histogram default parameters
      ParametersDescription hist_desc;
      // x-axis attributes
      hist_desc.add<std::vector<double> >("xbins", {0., 1.}).setDescription("x-axis bins definition");
      hist_desc.add<int>("nbins", 25).setDescription("Bins multiplicity for x-axis");
      hist_desc.add<int>("nbinsX", -1).setDescription("Bins multiplicity for x-axis");
      hist_desc.add<Limits>("xrange", Limits{}).setDescription("Minimum-maximum range for x-axis");
      // y-axis attributes
      hist_desc.add<std::vector<double> >("ybins", {0., 1.}).setDescription("y-axis bins definition");
      hist_desc.add<int>("nbinsY", 50).setDescription("Bins multiplicity for y-axis");
      hist_desc.add<Limits>("yrange", Limits{0., 1.}).setDescription("Minimum-maximum range for y-axis");
      hist_desc.add<bool>("log", false).setDescription("Plot logarithmic axis?");
      desc.addParametersDescriptionVector("histVariables", hist_desc)
          .setDescription("Histogram definition for 1/2 variable(s)");
      return desc;
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("text", TextHandler)
