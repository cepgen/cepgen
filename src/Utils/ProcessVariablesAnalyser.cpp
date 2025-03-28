/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/ProcessVariablesAnalyser.h"

using namespace cepgen;
using namespace cepgen::utils;

ProcessVariablesAnalyser::ProcessVariablesAnalyser(const proc::Process& proc, const ParametersList& params)
    : SteeredObject(params), proc_(proc) {
  for (const auto& var : proc_.mapped_variables_)
    if (auto hist = steer<ParametersList>(var.name); !hist.empty())
      hists_.insert(std::make_pair(var.name, new Hist1D(hist.set("name", var.name))));
    else
      hists_.insert(std::make_pair(
          var.name, new Hist1D(ParametersList().set("name", var.name).set("nbinsX", 50).set("xrange", var.limits))));
}

void ProcessVariablesAnalyser::feed(double weight) const {
  for (const auto& var : proc_.mapped_variables_)
    if (hists_.count(var.name))
      hists_.at(var.name)->fill(var.value, weight);
}

void ProcessVariablesAnalyser::analyse() const {
  auto drawer = DrawerFactory::get().build(steer<ParametersList>("drawer"));
  for (const auto& var : hists_)
    (void)drawer->draw(*var.second);
}

ParametersDescription ProcessVariablesAnalyser::description() {
  auto desc = ParametersDescription();
  ParametersDescription hist_desc;
  hist_desc.add("xbins", std::vector<double>{}).setDescription("x-axis bins definition");
  hist_desc.add("nbinsX", 25).setDescription("Bins multiplicity for x-axis");
  hist_desc.add("xrange", Limits{0., 1.}).setDescription("Minimum-maximum range for x-axis");
  desc.addParametersDescriptionVector("histVariables", hist_desc, {}).setDescription("Histogram definition");
  desc.add("drawer", DrawerFactory::get().describeParameters("root"));
  return desc;
}
