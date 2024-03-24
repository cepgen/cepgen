/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventBrowser.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::Generator mg;

  vector<string> vars;
  string input_card, plotter;
  int num_events;
  bool draw_grid, log;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("input,i", "input card", &input_card)
      .addArgument("vars", "variables to plot", &vars)
      .addOptionalArgument("num-events,n", "number of events to generate", &num_events, 100)
      .addOptionalArgument("draw-grid,g", "draw the x/y grid", &draw_grid, false)
      .addOptionalArgument("log,l", "logarithmic axis", &log, false)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "")
      .parse();

  mg.parseRunParameters(input_card);
  mg.runParameters().clearEventExportersSequence();

  // book all histograms
  map<string, unique_ptr<cepgen::utils::Hist1D> > h_var_hist;
  for (const auto& var : vars) {
    cepgen::Limits lim{0., 250.};
    if (var == "eta")
      lim = {-3., 3.};
    const auto parsed_var = cepgen::utils::replaceAll(var, ":", ",");
    h_var_hist[parsed_var].reset(new cepgen::utils::Hist1D(100, lim, cepgen::utils::sanitise(var)));
    h_var_hist[parsed_var]->xAxis().setLabel(parsed_var);
    h_var_hist[parsed_var]->yAxis().setLabel("d$\\sigma$/d" + parsed_var);
  }
  CG_DEBUG("main") << "Variables to be plotted: " << vars << ".";

  CG_LOG << "Process name: " << mg.runParameters().processName() << ".";

  cepgen::utils::EventBrowser browser;

  // generate the events and feed to the histogram(s)
  mg.generate(num_events, [&browser, &h_var_hist](const cepgen::Event& ev, unsigned long) {
    for (auto& var : h_var_hist)
      var.second->fill(browser.get(ev, var.first));
  });

  // normalise to cross section
  for (auto& var : h_var_hist)
    var.second->normalise(mg.crossSection() / num_events);

  // if a plotter is specified, draw histograms
  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    cepgen::utils::Drawer::Mode dm;
    if (log)
      dm |= cepgen::utils::Drawer::Mode::logy;
    if (draw_grid)
      dm |= cepgen::utils::Drawer::Mode::grid;

    for (const auto& var : h_var_hist)
      plt->draw(*var.second, dm);
  }

  return 0;
}
