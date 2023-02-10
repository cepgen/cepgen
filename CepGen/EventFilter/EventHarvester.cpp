/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/EventFilter/EventHarvester.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  EventHarvester::EventHarvester(const ParametersList& params)
      : EventExporter(params),
        browser_(new utils::EventBrowser),
        show_hists_(steer<bool>("show")),
        save_hists_(steer<bool>("save")),
        filename_(steer<std::string>("filename")) {
    // build the plotter object if specified
    const auto& plotter = steer<std::string>("plotter");
    if (!plotter.empty())
      drawer_ = DrawerFactory::get().build(plotter, params);

    // extract list of variables to be plotted in histogram
    const auto& hist_vars = steer<ParametersList>("histVariables");
    for (const auto& key : hist_vars.keys()) {
      const auto& vars = utils::split(key, ':');
      if (vars.size() < 1 || vars.size() > 2)
        throw CG_FATAL("EventHarvester") << "Invalid number of variables to correlate for '" << key << "'!";

      const auto& hvar = hist_vars.get<ParametersList>(key);
      const auto& log = hvar.get<bool>("log");
      auto name = utils::sanitise(key);
      if (vars.size() == 1) {  // 1D histogram
        const auto& xbins = hvar.get<std::vector<double> >("xbins");
        if (xbins.size() > 1)
          hists_.emplace_back(Hist1DInfo{vars.at(0), utils::Hist1D(xbins, name), log});
        else if (hvar.get<Limits>("xrange").valid()) {
          const auto& nbins = (hvar.get<int>("nbins") > 0 ? hvar.get<int>("nbins") : hvar.get<int>("nbinsX"));
          hists_.emplace_back(Hist1DInfo{vars.at(0), utils::Hist1D(nbins, hvar.get<Limits>("xrange"), name), log});
        } else {
          CG_WARNING("EventHarvester") << "Neither xrange nor xbins found in parameters for 1D plot of variable \""
                                       << vars.at(0) << "\".";
          continue;
        }
        auto& hist = hists_.rbegin()->hist;
        hist.xAxis().setLabel(vars.at(0));
        hist.yAxis().setLabel("d$\\sigma$/d" + vars.at(0) + " (pb/bin)");
      } else if (vars.size() == 2) {  // 2D histogram
        const auto &xbins = hvar.get<std::vector<double> >("xbins"), &ybins = hvar.get<std::vector<double> >("ybins");
        name = utils::sanitise(name);
        if (xbins.size() > 1 && ybins.size() > 1)
          hists2d_.emplace_back(Hist2DInfo{vars.at(0), vars.at(1), utils::Hist2D(xbins, ybins, name), log});
        else if (hvar.get<Limits>("xrange").valid()) {
          const auto& nbinsx = (hvar.get<int>("nbins") > 0 ? hvar.get<int>("nbins") : hvar.get<int>("nbinsX"));
          hists2d_.emplace_back(Hist2DInfo{
              vars.at(0),
              vars.at(1),
              utils::Hist2D(
                  nbinsx, hvar.get<Limits>("xrange"), hvar.get<int>("nbinsY"), hvar.get<Limits>("yrange"), name),
              log});
        } else {
          CG_WARNING("EventHarvester")
              << "Neither (x/y)range nor (x/y)bins found in parameters for 1D plot of variables \"" << vars << "\".";
          continue;
        }
        auto& hist = hists2d_.rbegin()->hist;
        hist.xAxis().setLabel(vars.at(0));
        hist.yAxis().setLabel(vars.at(1));
        hist.zAxis().setLabel("d$^2$$\\sigma$/d" + vars.at(0) + "/d" + vars.at(1) + " (pb/bin)");
      }
    }
    if (save_hists_ && !hists_.empty())
      file_.open(filename_);
  }

  EventHarvester::~EventHarvester() {
    //--- histograms printout
    if (!show_hists_ && !save_hists_)
      return;
    for (auto& h_var : hists_) {
      h_var.hist.scale(cross_section_ / (num_evts_ + 1));
      h_var.hist.setTitle(proc_name_);
      std::ostringstream os;
      if (drawer_)
        drawer_->draw(h_var.hist, h_var.log ? utils::Drawer::Mode::logy : utils::Drawer::Mode::none);
      if (show_hists_)
        CG_INFO("EventHarvester") << os.str();
      if (save_hists_)
        file_ << "\n" << os.str() << "\n";
    }
    for (auto& h_var : hists2d_) {
      std::ostringstream os;
      h_var.hist.setTitle(proc_name_);
      if (drawer_)
        drawer_->draw(h_var.hist,
                      utils::Drawer::Mode::grid | (h_var.log ? utils::Drawer::Mode::logz : utils::Drawer::Mode::none));
      if (show_hists_)
        CG_INFO("EventHarvester") << os.str();
      if (save_hists_)
        file_ << "\n" << os.str() << "\n";
    }
    if (save_hists_)
      CG_INFO("EventHarvester") << "Saved " << utils::s("histogram", hists_.size(), true) << " into \"" << filename_
                                << "\".";
  }

  void EventHarvester::initialise() {
    sqrts_ = runParameters().kinematics().incomingBeams().sqrtS();
    num_evts_ = 0ul;
    proc_name_ = ProcessFactory::get().describe(runParameters().processName());
    proc_name_ += ", \\sqrt{s} = " + utils::format("%g", sqrts_ * 1.e-3) + " TeV";
    if (save_hists_ && !hists_.empty())
      file_ << banner("#") << "\n";
  }

  void EventHarvester::operator<<(const Event& ev) {
    //--- increment the corresponding histograms
    for (auto& h_var : hists_)
      h_var.hist.fill(browser_->get(ev, h_var.var));
    for (auto& h_var : hists2d_)
      h_var.hist.fill(browser_->get(ev, h_var.var1), browser_->get(ev, h_var.var2));
    ++num_evts_;
  }

  ParametersDescription EventHarvester::description() {
    auto desc = EventExporter::description();
    desc.setDescription("Text-based histogramming tool");
    desc.add<std::string>("plotter", "").setDescription("Plotting algorithm to use");
    desc.add<std::string>("filename", "output.hists.txt").setDescription("Output filename for histogram dump");
    desc.add<bool>("show", true).setDescription("Show the histogram(s) at the end of the run?");
    desc.add<bool>("save", false).setDescription("Save the histogram(s) at the end of the run?");
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
}  // namespace cepgen
