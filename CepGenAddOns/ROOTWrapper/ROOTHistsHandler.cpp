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

namespace cepgen {
  /**
     * \brief Handler for the generic ROOT file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
  class ROOTHistsHandler final : public EventExporter {
  public:
    explicit ROOTHistsHandler(const ParametersList&);
    ~ROOTHistsHandler();

    static ParametersDescription description();

    void initialise() override {}
    void setCrossSection(const Value& cross_section) override { cross_section_ = cross_section; }
    void operator<<(const Event&) override;

  private:
    TFile file_;
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
        file_(steer<std::string>("filename").c_str(), "recreate"),
        variables_(steer<ParametersList>("variables")) {
    //--- extract list of variables/correlations to be plotted in histograms
    for (const auto& key : variables_.keys()) {
      const auto& vars = utils::split(key, ':');
      if (vars.size() < 1 || vars.size() > 3)
        throw CG_FATAL("ROOTHistsHandler") << "Invalid number of variables to correlate for '" << key << "'!";

      const auto& hvars = variables_.get<ParametersList>(key);
      int nbins_x = hvars.get<int>("nbinsX");
      if (hvars.get<int>("nbins") > 0)
        nbins_x = hvars.get<int>("nbins");
      const auto& xrange = hvars.get<Limits>("xrange");
      const bool profile = hvars.get<bool>("profile");
      if (vars.size() == 1) {  // 1D histogram
        auto title = hvars.get<std::string>("title");
        if (title.empty())
          title = utils::format("%s;%s;d#sigma/d(%s) (pb/bin)", key.c_str(), key.c_str(), key.c_str());
        hists1d_.emplace_back(
            std::make_pair(key, new TH1D(key.c_str(), title.c_str(), nbins_x, xrange.min(), xrange.max())));
        CG_INFO("ROOTHistsHandler") << "Booking a 1D histogram with " << utils::s("bin", nbins_x) << " in range "
                                    << xrange << " for \"" << key << "\".";
        continue;
      }
      const int nbins_y = hvars.get<int>("nbinsY");
      const auto& yrange = hvars.get<Limits>("yrange");
      if (vars.size() == 2) {  // 2D histogram / 1D profile
        auto title = hvars.get<std::string>("title");
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
              std::make_pair(vars, new TProfile(key.c_str(), title.c_str(), nbins_x, xrange.min(), xrange.max())));
          CG_INFO("ROOTHistsHandler") << "Booking a 1D profile with " << utils::s("bin", nbins_x, true) << " in range "
                                      << xrange << " for \"" << utils::merge(vars, " / ") << "\".";
        } else {
          hists2d_.emplace_back(std::make_pair(
              vars,
              new TH2D(
                  key.c_str(), title.c_str(), nbins_x, xrange.min(), xrange.max(), nbins_y, yrange.min(), yrange.max())));
          CG_INFO("ROOTHistsHandler") << "Booking a 2D correlation plot with "
                                      << utils::s("bin", nbins_x + nbins_y, true) << " in range x=" << xrange
                                      << " and y=" << yrange << " for \"" << utils::merge(vars, " / ") << "\".";
        }
        continue;
      }
      const int nbins_z = hvars.get<int>("nbinsZ");
      const auto& zrange = hvars.get<Limits>("zrange");
      if (vars.size() == 3) {  // 3D histogram
        auto title = hvars.get<std::string>("title");
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
          profiles2d_.emplace_back(std::make_pair(
              vars,
              new TProfile2D(
                  key.c_str(), title.c_str(), nbins_x, xrange.min(), xrange.max(), nbins_y, yrange.min(), yrange.max())));
          CG_INFO("ROOTHistsHandler") << "Booking a 2D profile with " << utils::s("bin", nbins_x + nbins_y, true)
                                      << " in range x=" << xrange << " and y=" << yrange << " for \""
                                      << utils::merge(vars, " / ") << "\".";
        } else {
          hists3d_.emplace_back(std::make_pair(vars,
                                               new TH3D(key.c_str(),
                                                        title.c_str(),
                                                        nbins_x,
                                                        xrange.min(),
                                                        xrange.max(),
                                                        nbins_y,
                                                        yrange.min(),
                                                        yrange.max(),
                                                        nbins_z,
                                                        zrange.min(),
                                                        zrange.max())));
          CG_INFO("ROOTHistsHandler") << "Booking a 3D correlation plot with "
                                      << utils::s("bin", nbins_x + nbins_y + nbins_z, true) << " in range x=" << xrange
                                      << ", y=" << yrange << ", and z=" << zrange << " for \""
                                      << utils::merge(vars, " / ") << "\".";
        }
        continue;
      }
    }
  }

  ROOTHistsHandler::~ROOTHistsHandler() {
    //--- finalisation of the output file
    for (const auto& hist : hists1d_)
      hist.second->Write(hist.first.c_str());
    for (const auto& hist : hists2d_)
      hist.second->Write(utils::merge(hist.first, "_vs_").c_str());
    for (const auto& hist : hists3d_)
      hist.second->Write(utils::merge(hist.first, "_vs_").c_str());
    for (const auto& hist : profiles1d_)
      hist.second->Write(utils::merge(hist.first, "_vs_").c_str());
    for (const auto& hist : profiles2d_)
      hist.second->Write(utils::merge(hist.first, "_vs_").c_str());
    // ROOT and its sumptuous memory management disallows the "delete" here
    file_.Close();
  }

  void ROOTHistsHandler::operator<<(const Event& ev) {
    //--- increment the corresponding histograms
    for (const auto& h_var : hists1d_)
      h_var.second->Fill(browser_.get(ev, h_var.first), (double)cross_section_);
    for (const auto& h_var : hists2d_)
      h_var.second->Fill(browser_.get(ev, h_var.first[0]), browser_.get(ev, h_var.first[1]), (double)cross_section_);
    for (const auto& h_var : hists3d_)
      h_var.second->Fill(browser_.get(ev, h_var.first[0]),
                         browser_.get(ev, h_var.first[1]),
                         browser_.get(ev, h_var.first[2]),
                         (double)cross_section_);
    for (const auto& h_var : profiles1d_)
      h_var.second->Fill(browser_.get(ev, h_var.first[0]), browser_.get(ev, h_var.first[1]), (double)cross_section_);
    for (const auto& h_var : profiles2d_)
      h_var.second->Fill(browser_.get(ev, h_var.first[0]),
                         browser_.get(ev, h_var.first[1]),
                         browser_.get(ev, h_var.first[2]),
                         (double)cross_section_);
  }

  ParametersDescription ROOTHistsHandler::description() {
    auto desc = EventExporter::description();
    desc.setDescription("ROOT histograming/profiling module");
    desc.add<std::string>("filename", "output.root").setDescription("Output filename");
    auto var_desc = ParametersDescription();
    var_desc.add<std::string>("title", "").setDescription("Variable description");
    var_desc.add<int>("nbins", -1);
    var_desc.add<int>("nbinsX", 10).setDescription("Bins multiplicity for x-axis");
    var_desc.add<Limits>("xrange", Limits{0., 1.}).setDescription("Minimum-maximum range for x-axis");
    var_desc.add<int>("nbinsY", 10).setDescription("Bins multiplicity for y-axis");
    var_desc.add<Limits>("yrange", Limits{0., 1.}).setDescription("Minimum-maximum range for y-axis");
    var_desc.add<int>("nbinsZ", 10).setDescription("Bins multiplicity for z-axis");
    var_desc.add<Limits>("zrange", Limits{0., 1.}).setDescription("Minimum-maximum range for z-axis");
    var_desc.add<bool>("profile", false);
    desc.addParametersDescriptionVector("variables", var_desc);
    return desc;
  }
}  // namespace cepgen

REGISTER_EXPORTER("root_hist", ROOTHistsHandler);
