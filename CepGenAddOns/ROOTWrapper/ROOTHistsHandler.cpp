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

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the generic ROOT file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class ROOTHistsHandler : public ExportModule {
    public:
      explicit ROOTHistsHandler(const ParametersList&);
      ~ROOTHistsHandler();
      static std::string description() { return "ROOT histograming/profiling module"; }

      void initialise(const Parameters&) override {}
      void setCrossSection(double cross_section, double) override { cross_section_ = cross_section; }
      void operator<<(const Event&) override;

    private:
      TFile file_;
      std::vector<std::pair<std::string, TH1*> > hists1d_;
      std::vector<std::pair<std::vector<std::string>, TH2*> > hists2d_;
      std::vector<std::pair<std::vector<std::string>, TH3*> > hists3d_;
      std::vector<std::pair<std::vector<std::string>, TProfile*> > profiles1d_;
      std::vector<std::pair<std::vector<std::string>, TProfile2D*> > profiles2d_;

      const ParametersList variables_;

      double cross_section_{1.};
      const utils::EventBrowser browser_;
    };

    ROOTHistsHandler::ROOTHistsHandler(const ParametersList& params)
        : ExportModule(params),
          file_(params.get<std::string>("filename", "output.root").c_str(), "recreate"),
          variables_(params.get<ParametersList>("variables")) {
      //--- extract list of variables/correlations to be plotted in histograms
      for (const auto& key : variables_.keys()) {
        const auto& vars = utils::split(key, ':');
        if (vars.size() < 1 || vars.size() > 3)
          throw CG_FATAL("ROOTHistsHandler") << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvars = variables_.get<ParametersList>(key);
        int nbins_x = hvars.get<int>("nbinsX", 10);
        nbins_x = hvars.get<int>("nbins", nbins_x);
        const auto& xrange = hvars.get<Limits>("xrange", Limits{0., 1.});
        const bool profile = hvars.get<bool>("profile", false);
        if (vars.size() == 1) {  // 1D histogram
          const auto title = hvars.get<std::string>(
              "title", utils::format("%s;%s;d#sigma/d(%s) (pb/bin)", key.c_str(), key.c_str(), key.c_str()));
          hists1d_.emplace_back(
              std::make_pair(key, new TH1D(key.c_str(), title.c_str(), nbins_x, xrange.min(), xrange.max())));
          CG_INFO("ROOTHistsHandler") << "Booking a 1D histogram with " << utils::s("bin", nbins_x) << " in range "
                                      << xrange << " for \"" << key << "\".";
          continue;
        }
        const int nbins_y = hvars.get<int>("nbinsY", 10);
        const auto& yrange = hvars.get<Limits>("yrange", Limits{0., 1.});
        if (vars.size() == 2) {  // 2D histogram / 1D profile
          const auto title =
              hvars.get<std::string>("title",
                                     utils::format("(%s / %s) correlation;%s;%s;d^{2}#sigma/d(%s)/d(%s) (pb/bin)",
                                                   vars[0].c_str(),
                                                   vars[1].c_str(),
                                                   vars[0].c_str(),
                                                   vars[1].c_str(),
                                                   vars[0].c_str(),
                                                   vars[1].c_str()));
          if (profile) {
            profiles1d_.emplace_back(
                std::make_pair(vars, new TProfile(key.c_str(), title.c_str(), nbins_x, xrange.min(), xrange.max())));
            CG_INFO("ROOTHistsHandler") << "Booking a 1D profile with " << utils::s("bin", nbins_x, true)
                                        << " in range " << xrange << " for \"" << utils::merge(vars, " / ") << "\".";
          } else {
            hists2d_.emplace_back(std::make_pair(vars,
                                                 new TH2D(key.c_str(),
                                                          title.c_str(),
                                                          nbins_x,
                                                          xrange.min(),
                                                          xrange.max(),
                                                          nbins_y,
                                                          yrange.min(),
                                                          yrange.max())));
            CG_INFO("ROOTHistsHandler") << "Booking a " << (profile ? "1D profile" : "2D correlation plot") << " with "
                                        << utils::s("bin", nbins_x + nbins_y, true) << " in range x=" << xrange
                                        << " and y=" << yrange << " for \"" << utils::merge(vars, " / ") << "\".";
          }
          continue;
        }
        const int nbins_z = hvars.get<int>("nbinsZ", 10);
        const auto& zrange = hvars.get<Limits>("zrange", Limits{0., 1.});
        if (vars.size() == 3) {  // 3D histogram
          const auto title = hvars.get<std::string>(
              "title",
              utils::format("(%s / %s / %s) correlation;%s;%s;%s;d^{3}#sigma/d(%s)/d(%s)/d(%s) (pb/bin)",
                            vars[0].c_str(),
                            vars[1].c_str(),
                            vars[2].c_str(),
                            vars[0].c_str(),
                            vars[1].c_str(),
                            vars[2].c_str(),
                            vars[0].c_str(),
                            vars[1].c_str(),
                            vars[2].c_str()));
          if (profile) {
            profiles2d_.emplace_back(std::make_pair(vars,
                                                    new TProfile2D(key.c_str(),
                                                                   title.c_str(),
                                                                   nbins_x,
                                                                   xrange.min(),
                                                                   xrange.max(),
                                                                   nbins_y,
                                                                   yrange.min(),
                                                                   yrange.max())));
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
                                        << utils::s("bin", nbins_x + nbins_y + nbins_z, true)
                                        << " in range x=" << xrange << ", y=" << yrange << ", and z=" << zrange
                                        << " for \"" << utils::merge(vars, " / ") << "\".";
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
        h_var.second->Fill(browser_.get(ev, h_var.first), cross_section_);
      for (const auto& h_var : hists2d_)
        h_var.second->Fill(browser_.get(ev, h_var.first[0]), browser_.get(ev, h_var.first[1]), cross_section_);
      for (const auto& h_var : hists3d_)
        h_var.second->Fill(browser_.get(ev, h_var.first[0]),
                           browser_.get(ev, h_var.first[1]),
                           browser_.get(ev, h_var.first[2]),
                           cross_section_);
      for (const auto& h_var : profiles1d_)
        h_var.second->Fill(browser_.get(ev, h_var.first[0]), browser_.get(ev, h_var.first[1]), cross_section_);
      for (const auto& h_var : profiles2d_)
        h_var.second->Fill(browser_.get(ev, h_var.first[0]),
                           browser_.get(ev, h_var.first[1]),
                           browser_.get(ev, h_var.first[2]),
                           cross_section_);
    }
  }  // namespace io
}  // namespace cepgen

REGISTER_IO_MODULE("root_hist", ROOTHistsHandler)
