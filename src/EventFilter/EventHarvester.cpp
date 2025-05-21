/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/String.h"

using namespace cepgen;
using namespace std::string_literals;

/// Generic per-event information output handler
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Jul 2019
class EventHarvester : public EventExporter {
public:
  explicit EventHarvester(const ParametersList& params) : EventExporter(params), browser_(new utils::EventBrowser) {
    if (const auto& plotter = steer<std::string>("plotter"); !plotter.empty())  // build the plotter object if specified
      drawer_ = DrawerFactory::get().build(plotter, params);

    // extract list of variables to be plotted in histogram
    const auto& hist_vars = steer<ParametersList>("histVariables");
    for (const auto& key : hist_vars.keys()) {
      const auto& vars = utils::split(key, ':');
      if (vars.size() < 1 || vars.size() > 2)
        throw CG_FATAL("EventHarvester") << "Invalid number of variables to correlate for '" << key << "'!";

      auto hvar = hist_vars.get<ParametersList>(key);
      const auto& log = hvar.get<bool>("log");
      auto name = utils::sanitise(key);
      if (vars.size() == 1) {  // 1D histogram
        auto hist = utils::Hist1D(hvar.set("name", name));
        hist.xAxis().setLabel(vars.at(0));
        hist.yAxis().setLabel("d$\\sigma$/d" + vars.at(0) + " (pb/bin)");
        hists_.emplace_back(Hist1DInfo{vars.at(0), hist, log});
      } else if (vars.size() == 2) {  // 2D histogram
        auto hist = utils::Hist2D(hvar.set("name", utils::sanitise(name)));
        hist.xAxis().setLabel(vars.at(0));
        hist.yAxis().setLabel(vars.at(1));
        hist.zAxis().setLabel("d$^2$$\\sigma$/d" + vars.at(0) + "/d" + vars.at(1) + " (pb/bin)");
        hists2d_.emplace_back(Hist2DInfo{vars.at(0), vars.at(1), hist, log});
      }
    }
  }
  ~EventHarvester() override {
    try {  // histogram printout
      for (auto& h_var : hists_) {
        h_var.hist.scale(cross_section_ / (num_events_ + 1));
        h_var.hist.setTitle(proc_name_);
        std::ostringstream os;
        if (drawer_)
          (void)drawer_->draw(h_var.hist, h_var.log ? utils::Drawer::Mode::logy : utils::Drawer::Mode::none);
      }
      for (auto& h_var : hists2d_) {
        std::ostringstream os;
        h_var.hist.setTitle(proc_name_);
        if (drawer_)
          (void)drawer_->draw(
              h_var.hist,
              utils::Drawer::Mode::grid | (h_var.log ? utils::Drawer::Mode::logz : utils::Drawer::Mode::none));
      }
    } catch (const Exception& error) {
      CG_ERROR("EventHarvester") << "Failed to save the histograms harvested in this run. Error received: "
                                 << error.what();
    }
  }

  static ParametersDescription description() {
    auto desc = EventExporter::description();
    desc.setDescription("Event-based histogramming tool");
    desc.add("plotter"s, ""s).setDescription("Plotting algorithm to use");
    desc.addParametersDescriptionVector("histVariables"s, utils::Hist2D::description(), {})
        .setDescription("Histogram definition for 1/2 variable(s)");
    return desc;
  }

  inline void setCrossSection(const Value& cross_section) override { cross_section_ = cross_section; }
  bool operator<<(const Event& event) override {
    // increment the corresponding histograms
    for (auto& h_var : hists_)
      h_var.hist.fill(browser_->get(event, h_var.var));
    for (auto& h_var : hists2d_)
      h_var.hist.fill(browser_->get(event, h_var.var1), browser_->get(event, h_var.var2));
    ++num_events_;
    return true;
  }

private:
  void initialise() override {
    num_events_ = 0ul;
    proc_name_ = ProcessFactory::get().describe(runParameters().processName());
    proc_name_ +=
        ", \\sqrt{s} = " + utils::format("%g", runParameters().kinematics().incomingBeams().sqrtS() * 1.e-3) + " TeV";
  }

  const std::unique_ptr<utils::EventBrowser> browser_;  ///< Event string-to-quantity extraction tool
  std::unique_ptr<utils::Drawer> drawer_;               ///< Drawing utility
  Value cross_section_{1., 0.};                         ///< Cross-section value, in pb
  unsigned long num_events_{0ul};                       ///< Number of events processed
  std::string proc_name_;                               ///< Name of the physics process

  /// 1D histogram definition
  struct Hist1DInfo {
    std::string var;
    utils::Hist1D hist;
    bool log;
  };
  std::vector<Hist1DInfo> hists_;  ///< List of 1D histograms
  /// 2D histogram definition
  struct Hist2DInfo {
    std::string var1;
    std::string var2;
    utils::Hist2D hist;
    bool log;
  };
  std::vector<Hist2DInfo> hists2d_;  ///< List of 2D histograms
};
REGISTER_EXPORTER("eventHarvester", EventHarvester);

struct TextHarvester final : EventHarvester {
  using EventHarvester::EventHarvester;
  static ParametersDescription description() {
    auto desc = EventHarvester::description();
    desc.setDescription("Text-based event harvester");
    desc.add("plotter", "text"s);
    return desc;
  }
};
REGISTER_EXPORTER("text", TextHarvester);
