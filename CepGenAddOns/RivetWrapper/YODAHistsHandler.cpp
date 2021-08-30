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

#include <limits>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ExportModule.h"
#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Utils/Limits.h"
#include "CepGen/Utils/String.h"
#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/WriterAIDA.h"
#include "YODA/WriterFLAT.h"
#include "YODA/WriterYODA.h"

namespace cepgen {
  namespace io {
    /**
     * \brief Handler for the generic YODA file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     * \tparam T YODA writer type
     */
    template <typename T>
    class YODAHistsHandler : public ExportModule {
    public:
      explicit YODAHistsHandler(const ParametersList&);
      ~YODAHistsHandler();
      static std::string description() { return "YODA histograms/profiles file output module"; }

      void initialise(const Parameters&) override {}
      void setCrossSection(double cross_section, double) override { cross_section_ = cross_section; }
      void operator<<(const Event&) override;

    private:
      std::ofstream file_;
      std::vector<std::pair<std::string, YODA::Histo1D> > hists1d_;
      std::vector<std::pair<std::vector<std::string>, YODA::Histo2D> > hists2d_;
      std::vector<std::pair<std::vector<std::string>, YODA::Profile1D> > profiles1d_;
      std::vector<std::pair<std::vector<std::string>, YODA::Profile2D> > profiles2d_;
      YODA::Counter weight_cnt_;
      const ParametersList variables_;

      double cross_section_{1.};
      const utils::EventBrowser browser_;
    };

    template <typename T>
    YODAHistsHandler<T>::YODAHistsHandler(const ParametersList& params)
        : ExportModule(params),
          file_(params.get<std::string>("filename", "output.yoda")),
          variables_(params.get<ParametersList>("variables")) {
      //--- extract list of variables/correlations to be plotted in histograms
      for (const auto& key : variables_.keys()) {
        const auto& vars = utils::split(key, ':');
        if (vars.size() < 1 || vars.size() > 3)
          throw CG_FATAL("YODAHistsHandler") << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvars = variables_.get<ParametersList>(key);
        int nbins_x = hvars.get<int>("nbinsX", 10);
        nbins_x = hvars.get<int>("nbins", nbins_x);
        const auto& xrange = hvars.get<Limits>("xrange", Limits{0., 1.});
        const bool profile = hvars.get<bool>("profile", false);
        if (vars.size() == 1) {  // 1D histogram
          const auto title = utils::format("d(sigma)/d(%s) (pb/bin)", key.c_str());
          hists1d_.emplace_back(std::make_pair(key, YODA::Histo1D(nbins_x, xrange.min(), xrange.max(), key, title)));
          CG_INFO("YODAHistsHandler") << "Booking a histogram with " << utils::s("bin", nbins_x) << " in range "
                                      << xrange << " for \"" << vars[0] << "\".";
          continue;
        }
        const int nbins_y = hvars.get<int>("nbinsY", 10);
        const auto& yrange = hvars.get<Limits>("yrange", Limits{0., 1.});
        if (vars.size() == 2) {  // 2D histogram
          const auto title = utils::format("d^2(sigma)/d(%s)/d(%s) (pb/bin)", vars[0].c_str(), vars[1].c_str());
          if (profile) {
            profiles1d_.emplace_back(
                std::make_pair(vars, YODA::Profile1D(nbins_x, xrange.min(), xrange.max(), key, title)));
            CG_INFO("YODAHistsHandler") << "Booking a 1D profile with " << utils::s("bin", nbins_x)
                                        << " in range x=" << xrange << " for \"" << utils::merge(vars, " / ") << "\".";
          } else {
            hists2d_.emplace_back(std::make_pair(
                vars,
                YODA::Histo2D(nbins_x, xrange.min(), xrange.max(), nbins_y, yrange.min(), yrange.max(), key, title)));
            CG_INFO("YODAHistsHandler") << "Booking a 2D correlation plot with " << utils::s("bin", nbins_x + nbins_y)
                                        << " in range x=" << xrange << " and y=" << yrange << " for \""
                                        << utils::merge(vars, " / ") << "\".";
          }
          continue;
        }
        if (vars.size() == 3 && profile) {
          const auto title = utils::format("(%s / %s / %s) correlation;%s;%s;%s;d^{3}#sigma/d(%s)/d(%s)/d(%s) (pb/bin)",
                                           vars[0].c_str(),
                                           vars[1].c_str(),
                                           vars[2].c_str(),
                                           vars[0].c_str(),
                                           vars[1].c_str(),
                                           vars[2].c_str(),
                                           vars[0].c_str(),
                                           vars[1].c_str(),
                                           vars[2].c_str());
          profiles2d_.emplace_back(std::make_pair(
              vars,
              YODA::Profile2D(nbins_x, xrange.min(), xrange.max(), nbins_y, yrange.min(), yrange.max(), key, title)));
          CG_INFO("YODAHistsHandler") << "Booking a 2D profile"
                                      << " with " << utils::s("bin", nbins_x + nbins_y, true)
                                      << " in range x=" << xrange << " and y=" << yrange << " for \""
                                      << utils::merge(vars, " / ") << "\".";
          continue;
        }
      }
    }

    template <typename T>
    YODAHistsHandler<T>::~YODAHistsHandler() {
      std::vector<const YODA::AnalysisObject*> obj;
      //--- finalisation of the output file
      for (const auto& hist : hists1d_)
        obj.emplace_back(&hist.second);
      for (const auto& hist : hists2d_)
        obj.emplace_back(&hist.second);
      for (const auto& hist : profiles1d_)
        obj.emplace_back(&hist.second);
      for (const auto& hist : profiles2d_)
        obj.emplace_back(&hist.second);
      obj.emplace_back(&weight_cnt_);
      T::write(file_, obj);
    }

    template <typename T>
    void YODAHistsHandler<T>::operator<<(const Event& ev) {
      //--- increment the corresponding histograms
      for (auto& h_var : hists1d_)
        h_var.second.fillBin(browser_.get(ev, h_var.first), cross_section_);
      for (auto& h_var : hists2d_)
        h_var.second.fillBin(browser_.get(ev, h_var.first[0]), browser_.get(ev, h_var.first[1]), cross_section_);
      for (auto& h_var : profiles1d_)
        h_var.second.fill(browser_.get(ev, h_var.first[0]), browser_.get(ev, h_var.first[1]), cross_section_);
      for (auto& h_var : profiles2d_)
        h_var.second.fill(browser_.get(ev, h_var.first[0]),
                          browser_.get(ev, h_var.first[1]),
                          browser_.get(ev, h_var.first[2]),
                          cross_section_);
      weight_cnt_.fill(ev.weight);
    }
  }  // namespace io
}  // namespace cepgen

typedef cepgen::io::YODAHistsHandler<YODA::WriterYODA> YodaOutputHandler;
typedef cepgen::io::YODAHistsHandler<YODA::WriterAIDA> YodaAidaOutputHandler;
typedef cepgen::io::YODAHistsHandler<YODA::WriterFLAT> YodaFlatOutputHandler;
REGISTER_IO_MODULE("yoda", YodaOutputHandler)
REGISTER_IO_MODULE("yoda_aida", YodaAidaOutputHandler)
REGISTER_IO_MODULE("yoda_flat", YodaFlatOutputHandler)
