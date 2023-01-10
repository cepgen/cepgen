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

#include <fstream>

#include "CepGen/Core/EventExporter.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Utils/Histogram.h"

namespace cepgen {
  class Event;
  class Parameters;
  namespace utils {
    class Drawer;
  }
  /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
  class IntegratedEventVariablesHandler : public EventExporter {
  public:
    explicit IntegratedEventVariablesHandler(const ParametersList&);
    ~IntegratedEventVariablesHandler();

    static ParametersDescription description();

    void initialise() override;
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
    /// Name of the physics process
    std::string proc_name_;

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
}  // namespace cepgen
