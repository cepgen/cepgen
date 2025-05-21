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
  explicit EventHarvester(const ParametersList& params) : EventExporter(params) {
    if (const auto& plotter = steer<std::string>("plotter"); !plotter.empty())  // build the plotter object if specified
      drawer_ = DrawerFactory::get().build(plotter, params);

    const auto& hist_vars = steer<ParametersList>("histVariables");
    for (const auto& key : hist_vars.keys()) {  // extract list of variables to be plotted in histogram
      auto hvar = hist_vars.get<ParametersList>(key);
      const auto& log = hvar.get<bool>("log");
      auto name = utils::sanitise(key);
      if (const auto& vars = utils::split(key, ':'); vars.size() == 1) {  // 1D histogram
        auto hist = utils::Hist1D(hvar.set("name", name));
        hist.xAxis().setLabel(vars.at(0));
        hist.yAxis().setLabel("d$\\sigma$/d" + vars.at(0) + " (pb/bin)");
        hists1d_.emplace_back(Hist1DInfo{vars.at(0), hist, log});
      } else if (vars.size() == 2) {  // 2D histogram
        auto hist = utils::Hist2D(hvar.set("name", utils::sanitise(name)));
        hist.xAxis().setLabel(vars.at(0));
        hist.yAxis().setLabel(vars.at(1));
        hist.zAxis().setLabel("d${}^2\\sigma$/d" + vars.at(0) + "/d" + vars.at(1) + " (pb/bin)");
        hists2d_.emplace_back(Hist2DInfo{vars.at(0), vars.at(1), hist, log});
      } else
        throw CG_FATAL("EventHarvester") << "Invalid number of variables to correlate for '" << key << "'.";
    }
  }
  ~EventHarvester() override {
    try {  // histogram printout
      for (auto& info : hists1d_) {
        info.histogram.scale(cross_section_ / num_events_);
        info.histogram.setTitle(proc_name_);
        if (drawer_)
          (void)drawer_->draw(info.histogram, info.log_y ? utils::Drawer::Mode::logy : utils::Drawer::Mode::none);
      }
      for (auto& info : hists2d_) {
        info.histogram.setTitle(proc_name_);
        if (drawer_)
          (void)drawer_->draw(
              info.histogram,
              utils::Drawer::Mode::grid | (info.log_z ? utils::Drawer::Mode::logz : utils::Drawer::Mode::none));
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

  void setCrossSection(const Value& cross_section) override { cross_section_ = cross_section; }
  bool operator<<(const Event& event) override {
    // increment the corresponding histograms
    for (auto& info : hists1d_)
      info.histogram.fill(browser_.get(event, info.variable));
    for (auto& info : hists2d_)
      info.histogram.fill(browser_.get(event, info.variable1), browser_.get(event, info.variable2));
    ++num_events_;
    return true;
  }

private:
  void initialise() override {
    num_events_ = 0ul;
    proc_name_ = ProcessFactory::get().describe(runParameters().processName());
    proc_name_ +=
        ", $\\sqrt{s} =$ " + utils::format("%g TeV", runParameters().kinematics().incomingBeams().sqrtS() * 1.e-3);
  }

  const utils::EventBrowser browser_;      ///< Event string-to-quantity extraction tool
  std::unique_ptr<utils::Drawer> drawer_;  ///< Drawing utility

  Value cross_section_{1., 0.};    ///< Cross-section value, in pb
  unsigned long num_events_{0ul};  ///< Number of events processed
  std::string proc_name_;          ///< Name of the physics process
  struct Hist1DInfo {
    std::string variable;
    utils::Hist1D histogram;
    bool log_y;
  };  ///< 1D histogram definition
  std::vector<Hist1DInfo> hists1d_;  ///< List of 1D histograms
  struct Hist2DInfo {
    std::string variable1;
    std::string variable2;
    utils::Hist2D histogram;
    bool log_z;
  };  ///< 2D histogram definition
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
