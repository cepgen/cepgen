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
using namespace std::string_literals;

ProcessVariablesAnalyser::ProcessVariablesAnalyser(const proc::Process& proc, const ParametersList& params)
    : SteeredObject(params), proc_(proc), drawer_(DrawerFactory::get().build(steer<ParametersList>("drawer"s))) {
  for (const auto& var : proc_.mapped_variables_)
    if (auto hist = steer<ParametersList>(var.name); !hist.empty())
      hists_.insert(std::make_pair(var.name, new Hist1D(hist.set("name", var.name))));
    else
      hists_.insert(std::make_pair(
          var.name, new Hist1D(ParametersList().set("name"s, var.name).set("nbinsX"s, 50).set("xrange"s, var.limits))));
}

void ProcessVariablesAnalyser::feed(double weight) const {
  for (const auto& var : proc_.mapped_variables_)
    if (hists_.count(var.name))
      hists_.at(var.name)->fill(var.value, weight);
}

void ProcessVariablesAnalyser::analyse() const {
  for (const auto& [variable_name, histogram] : hists_)
    (void)drawer_->draw(*histogram);
}

ParametersDescription ProcessVariablesAnalyser::description() {
  auto desc = ParametersDescription();
  desc.addParametersDescriptionVector("histVariables", Hist1D::description(), {})
      .setDescription("Histogram definition");
  desc.add("drawer", DrawerFactory::get().describeParameters("root"));
  return desc;
}
