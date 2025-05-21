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

#ifndef CepGen_EventFilter_EventHarvester_h
#define CepGen_EventFilter_EventHarvester_h

#include <fstream>

#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Utils/Histogram.h"

namespace cepgen {
  class Event;
}
namespace cepgen::utils {
  class Drawer;
  class EventBrowser;
}  // namespace cepgen::utils

namespace cepgen {
  /// Generic text file output handler
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jul 2019
  class EventHarvester : public EventExporter {
  public:
    explicit EventHarvester(const ParametersList&);
    ~EventHarvester() override;

    static ParametersDescription description();

    inline void setCrossSection(const Value& cross_section) override { cross_section_ = cross_section; }
    bool operator<<(const Event&) override;

  private:
    void initialise() override;

    const std::unique_ptr<utils::EventBrowser> browser_;  ///< Event string-to-quantity extaction tool
    const bool show_hists_;                               ///< Display histograms after the run
    const bool save_hists_;                               ///< Save histograms into the output file after the run
    const std::string filename_;                          ///< Output file path

    std::ofstream file_;                     ///< Output file where all information is stored
    std::unique_ptr<utils::Drawer> drawer_;  ///< Drawing utility
    Value cross_section_{1., 0.};            ///< Cross-section value, in pb
    unsigned long num_events_{0ul};          ///< Number of events processed
    std::string proc_name_;                  ///< Name of the physics process

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
}  // namespace cepgen

#endif
