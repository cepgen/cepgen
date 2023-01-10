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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

/** Example executable for CepGen
 * - loads the steering card variables into the environment,
 * - launches the cross-section computation and the events generation (if requested).
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main(int argc, char* argv[]) {
  string input_card;
  int num_events;
  bool list_mods, safe_mode;
  vector<string> addons, outputs;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addOptionalArgument("config,i", "path to the configuration file", &input_card)
      .addOptionalArgument("num-events,n", "number of events to generate", &num_events, -1)
      .addOptionalArgument("list-modules,l", "list all runtime modules", &list_mods, false)
      .addOptionalArgument("add-ons,a", "external runtime plugin", &addons)
      .addOptionalArgument("output,o", "additional output module(s)", &outputs)
      .addOptionalArgument("safe-mode,s", "safe mode", &safe_mode, false)
      .parse();

  cepgen::Generator gen(safe_mode);

  //--- first start by defining the generator object
  for (const auto& lib : addons)
    try {
      cepgen::loadLibrary(lib);
    } catch (const cepgen::Exception& e) {
      e.dump();
    }

  //--- if modules listing is requested
  if (list_mods) {
    cepgen::dumpModules();
    return 0;
  }

  //--- no steering card nor additional flags found
  if (input_card.empty() && parser.extra_config().empty()) {
    CG_WARNING("main") << "Neither input card nor configuration word provided!\n\n " << parser.help_message();
    return 0;
  }
  //--- parse the steering card
  if (!input_card.empty())
    gen.setParameters(cepgen::card::Handler::parse(input_card));
  //--- parse the additional flags
  if (!parser.extra_config().empty())
    gen.setParameters(cepgen::card::CardsHandlerFactory::get()
                          .build(".cmd", cepgen::ParametersList().set<vector<string> >("args", parser.extra_config()))
                          ->parse(string(), gen.parametersPtr()));

  cepgen::utils::AbortHandler();

  try {
    auto& params = gen.parametersRef();
    if (num_events >= 0)  // user specified a number of events to generate
      params.generation().setMaxGen(num_events);

    if (params.generation().enabled() && !outputs.empty())
      for (const auto& output : outputs)
        gen.parametersRef().addEventExporter(cepgen::EventExporterFactory::get().build(output));

    //--- list all parameters
    CG_LOG << gen.parameters();

    //--- let there be a cross-section...
    double xsec = 0., err = 0.;
    gen.computeXsection(xsec, err);

    if (params.generation().enabled())
      //--- events generation starts here
      // (one may use a callback function)
      gen.generate();
  } catch (const cepgen::utils::RunAbortedException&) {
    CG_DEBUG("main") << "Run aborted!";
  } catch (const cepgen::Exception& e) {
    e.dump();
  } catch (const exception& e) {
    CG_FATAL("main") << "Other exception caught!\n\t" << e.what();
  }

  return 0;
}
