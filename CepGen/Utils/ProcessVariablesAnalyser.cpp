/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/ProcessVariablesAnalyser.h"

namespace cepgen {
  namespace utils {
    ProcessVariablesAnalyser::ProcessVariablesAnalyser(const ParametersList& params) : SteeredObject(params) {}

    ProcessVariablesAnalyser& ProcessVariablesAnalyser::get(const ParametersList& params) {
      static ProcessVariablesAnalyser analyser(params);
      return analyser;
    }

    void ProcessVariablesAnalyser::reset(const proc::Process& proc) {
      hists_.clear();
      for (const auto& var : proc.mapped_variables_) {
        if (const auto& hist = steer<ParametersList>(var.name); !hist.empty()) {
          if (const auto& xbins = hist.get<std::vector<double> >("xbins"); xbins.size() > 1)
            hists_.insert(std::make_pair(var.name, Hist1D(xbins, var.name)));
          else if (hist.get<Limits>("xrange").valid()) {
            const auto& nbins = (hist.get<int>("nbins") > 0 ? hist.get<int>("nbins") : hist.get<int>("nbinsX"));
            hists_.insert(std::make_pair(var.name, Hist1D(nbins, hist.get<Limits>("xrange"), var.name)));
          }
        }
      }
    }

    void ProcessVariablesAnalyser::analyse(proc::Process& proc, double weight) {
      for (const auto& var : proc.mapped_variables_) {
        if (hists_.count(var.name))
          hists_.at(var.name).fill(var.value, weight);
      }
    }

    ParametersDescription ProcessVariablesAnalyser::description() {
      auto desc = ParametersDescription();
      ParametersDescription hist_desc;
      hist_desc.add<std::vector<double> >("xbins", {}).setDescription("x-axis bins definition");
      hist_desc.add<int>("nbins", 25).setDescription("Bins multiplicity for x-axis");
      hist_desc.add<int>("nbinsX", -1).setDescription("Bins multiplicity for x-axis");
      hist_desc.add<Limits>("xrange", Limits{0., 1.}).setDescription("Minimum-maximum range for x-axis");
      desc.addParametersDescriptionVector("histVariables", hist_desc, {}).setDescription("Histogram definition");
      return desc;
    }
  }  // namespace utils
}  // namespace cepgen
