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

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Utils/Value.h"

using namespace cepgen;
using namespace std::string_literals;

/// Handler for the generic ROOT file output
/// \author Laurent Forthomme <laurent.forthomme@cern.ch>
/// \date Jul 2019
class ROOTHistsHandler final : public EventExporter {
public:
  explicit ROOTHistsHandler(const ParametersList&);
  ~ROOTHistsHandler() override {
    // finalisation of the output file
    for (const auto& [names, hist] : hists1d_)
      hist->Write(names.c_str());
    for (const auto& [names, hist] : hists2d_)
      hist->Write(utils::merge(names, "_vs_").c_str());
    for (const auto& [names, hist] : hists3d_)
      hist->Write(utils::merge(names, "_vs_").c_str());
    for (const auto& [names, hist] : profiles1d_)
      hist->Write(utils::merge(names, "_vs_").c_str());
    for (const auto& [names, hist] : profiles2d_)
      hist->Write(utils::merge(names, "_vs_").c_str());
    file_->Close();  // ROOT and its sumptuous memory management disallow the "delete" here
  }

  static ParametersDescription description() {
    auto desc = EventExporter::description();
    desc.setDescription("ROOT histogramming/profiling module");
    desc.add("filename"s, "output.root"s).setDescription("Output filename");
    auto var_desc = ParametersDescription();
    var_desc.add("title"s, ""s).setDescription("Variable description");
    var_desc.add("nbins"s, -1);
    var_desc.add("nbinsX"s, 10).setDescription("Bins multiplicity for x-axis");
    var_desc.add("xrange"s, Limits{0., 1.}).setDescription("Minimum-maximum range for x-axis");
    var_desc.add("nbinsY"s, 10).setDescription("Bins multiplicity for y-axis");
    var_desc.add("yrange"s, Limits{0., 1.}).setDescription("Minimum-maximum range for y-axis");
    var_desc.add("nbinsZ"s, 10).setDescription("Bins multiplicity for z-axis");
    var_desc.add("zrange"s, Limits{0., 1.}).setDescription("Minimum-maximum range for z-axis");
    var_desc.add("profile"s, false);
    desc.addParametersDescriptionVector("variables", var_desc);
    return desc;
  }

  void setCrossSection(const Value& cross_section) override { cross_section_ = cross_section; }
  bool operator<<(const Event& event) override {
    // increment the corresponding histograms
    for (const auto& h_var : hists1d_)
      h_var.second->Fill(browser_.get(event, h_var.first), cross_section_);
    for (const auto& h_var : hists2d_)
      h_var.second->Fill(browser_.get(event, h_var.first[0]), browser_.get(event, h_var.first[1]), cross_section_);
    for (const auto& h_var : hists3d_)
      h_var.second->Fill(browser_.get(event, h_var.first[0]),
                         browser_.get(event, h_var.first[1]),
                         browser_.get(event, h_var.first[2]),
                         cross_section_);
    for (const auto& h_var : profiles1d_)
      h_var.second->Fill(browser_.get(event, h_var.first[0]), browser_.get(event, h_var.first[1]), cross_section_);
    for (const auto& h_var : profiles2d_)
      h_var.second->Fill(browser_.get(event, h_var.first[0]),
                         browser_.get(event, h_var.first[1]),
                         browser_.get(event, h_var.first[2]),
                         cross_section_);
    return true;
  }

private:
  const std::unique_ptr<TFile> file_;
  std::vector<std::pair<std::string, TH1*> > hists1d_;
  std::vector<std::pair<std::vector<std::string>, TH2*> > hists2d_;
  std::vector<std::pair<std::vector<std::string>, TH3*> > hists3d_;
  std::vector<std::pair<std::vector<std::string>, TProfile*> > profiles1d_;
  std::vector<std::pair<std::vector<std::string>, TProfile2D*> > profiles2d_;

  const ParametersList variables_;

  Value cross_section_{1., 0.};
  const utils::EventBrowser browser_;
};

ROOTHistsHandler::ROOTHistsHandler(const ParametersList& params)
    : EventExporter(params),
      file_(TFile::Open(steer<std::string>("filename"s).c_str(), "recreate")),
      variables_(steer<ParametersList>("variables"s)) {
  //--- extract list of variables/correlations to be plotted in histograms
  for (const auto& key : variables_.keys()) {
    const auto& vars = utils::split(key, ':');
    if (vars.size() < 1 || vars.size() > 3)
      throw CG_FATAL("ROOTHistsHandler") << "Invalid number of variables to correlate for '" << key << "'!";

    const auto& variable = variables_.get<ParametersList>(key);
    auto num_bins_x = variable.get<int>("nbinsX"s);
    if (variable.get<int>("nbins"s) > 0)
      num_bins_x = variable.get<int>("nbins"s);
    const auto& x_range = variable.get<Limits>("xrange"s);
    const bool profile = variable.get<bool>("profile"s);
    if (vars.size() == 1) {  // 1D histogram
      auto title = variable.get<std::string>("title"s);
      if (title.empty())
        title = utils::format("%s;%s;d#sigma/d(%s) (pb/bin)", key.c_str(), key.c_str(), key.c_str());
      hists1d_.emplace_back(
          std::make_pair(key, new TH1D(key.c_str(), title.c_str(), num_bins_x, x_range.min(), x_range.max())));
      CG_INFO("ROOTHistsHandler") << "Booking a 1D histogram with " << utils::s("bin", num_bins_x) << " in range "
                                  << x_range << " for \"" << key << "\".";
      continue;
    }
    const auto num_bins_y = variable.get<int>("nbinsY"s);
    const auto& y_range = variable.get<Limits>("yrange"s);
    if (vars.size() == 2) {  // 2D histogram / 1D profile
      auto title = variable.get<std::string>("title"s);
      if (title.empty())
        title = utils::format("(%s / %s) correlation;%s;%s;d^{2}#sigma/d(%s)/d(%s) (pb/bin)",
                              vars[0].c_str(),
                              vars[1].c_str(),
                              vars[0].c_str(),
                              vars[1].c_str(),
                              vars[0].c_str(),
                              vars[1].c_str());
      if (profile) {
        profiles1d_.emplace_back(
            std::make_pair(vars, new TProfile(key.c_str(), title.c_str(), num_bins_x, x_range.min(), x_range.max())));
        CG_INFO("ROOTHistsHandler") << "Booking a 1D profile with " << utils::s("bin", num_bins_x, true) << " in range "
                                    << x_range << " for \"" << utils::merge(vars, " / ") << "\".";
      } else {
        hists2d_.emplace_back(std::make_pair(vars,
                                             new TH2D(key.c_str(),
                                                      title.c_str(),
                                                      num_bins_x,
                                                      x_range.min(),
                                                      x_range.max(),
                                                      num_bins_y,
                                                      y_range.min(),
                                                      y_range.max())));
        CG_INFO("ROOTHistsHandler") << "Booking a 2D correlation plot with "
                                    << utils::s("bin", num_bins_x + num_bins_y, true) << " in range x=" << x_range
                                    << " and y=" << y_range << " for \"" << utils::merge(vars, " / ") << "\".";
      }
      continue;
    }
    const int num_bins_z = variable.get<int>("nbinsZ"s);
    const auto& z_range = variable.get<Limits>("zrange"s);
    if (vars.size() == 3) {  // 3D histogram
      auto title = variable.get<std::string>("title"s);
      if (title.empty())
        title = utils::format("(%s / %s / %s) correlation;%s;%s;%s;d^{3}#sigma/d(%s)/d(%s)/d(%s) (pb/bin)",
                              vars[0].c_str(),
                              vars[1].c_str(),
                              vars[2].c_str(),
                              vars[0].c_str(),
                              vars[1].c_str(),
                              vars[2].c_str(),
                              vars[0].c_str(),
                              vars[1].c_str(),
                              vars[2].c_str());
      if (profile) {
        profiles2d_.emplace_back(std::make_pair(vars,
                                                new TProfile2D(key.c_str(),
                                                               title.c_str(),
                                                               num_bins_x,
                                                               x_range.min(),
                                                               x_range.max(),
                                                               num_bins_y,
                                                               y_range.min(),
                                                               y_range.max())));
        CG_INFO("ROOTHistsHandler") << "Booking a 2D profile with " << utils::s("bin", num_bins_x + num_bins_y, true)
                                    << " in range x=" << x_range << " and y=" << y_range << " for \""
                                    << utils::merge(vars, " / ") << "\".";
      } else {
        hists3d_.emplace_back(std::make_pair(vars,
                                             new TH3D(key.c_str(),
                                                      title.c_str(),
                                                      num_bins_x,
                                                      x_range.min(),
                                                      x_range.max(),
                                                      num_bins_y,
                                                      y_range.min(),
                                                      y_range.max(),
                                                      num_bins_z,
                                                      z_range.min(),
                                                      z_range.max())));
        CG_INFO("ROOTHistsHandler") << "Booking a 3D correlation plot with "
                                    << utils::s("bin", num_bins_x + num_bins_y + num_bins_z, true)
                                    << " in range x=" << x_range << ", y=" << y_range << ", and z=" << z_range
                                    << " for \"" << utils::merge(vars, " / ") << "\".";
      }
    }
  }
}
REGISTER_EXPORTER("root_hist", ROOTHistsHandler);
