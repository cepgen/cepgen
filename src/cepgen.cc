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

#include "CepGen/Cards/Handler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/CardsHandlerFactory.h"
#include "CepGen/Modules/DocumentationGeneratorFactory.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/DocumentationGenerator.h"

using namespace std;

/** Example executable for CepGen
 * - loads the steering card variables into the environment,
 * - launches the cross-section computation and the events generation (if requested).
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 */
int main(int argc, char* argv[]) {
  string input_card;
  int num_events;
  bool list_mods;
  vector<string> outputs;

  cepgen::Generator gen;  // first start by defining the generator object
  cepgen::ArgumentsParser parser(argc, argv);
  parser.addOptionalArgument("config,i", "path to the configuration file", &input_card)
      .addOptionalArgument("num-events,n", "number of events to generate", &num_events, -1)
      .addOptionalArgument("list-modules,l", "list all runtime modules", &list_mods, false)
      .addOptionalArgument("output,o", "additional output module(s)", &outputs)
      .parse();

  if (list_mods) {  // modules listing is requested ; dump and exit
    auto doc_dump =
        cepgen::DocumentationGeneratorFactory::get().build("text", cepgen::ParametersList().set("light", true));
    CG_LOG << doc_dump->describe();
    return EXIT_SUCCESS;
  }

  if (input_card.empty() && parser.extra_config().empty()) {  // no steering card nor additional flags found
    CG_ERROR("main") << "Neither input card nor configuration word provided!\n\n " << parser.help_message();
    return EXIT_FAILURE;
  }
  if (!input_card.empty())
    gen.parseRunParameters(input_card);  // parse the steering card
  if (!parser.extra_config().empty()) {  // parse any additional flags
    auto args_handler = cepgen::CardsHandlerFactory::get().build(".cmd");
    args_handler->setRunParameters(&gen.runParameters());
    args_handler->parseCommands(parser.extra_config());
    gen.setRunParameters(args_handler->runParameters());
  }

  cepgen::utils::AbortHandler();  // allow C-c to be triggered during a run

  try {
    auto& params = gen.runParameters();
    if (num_events >= 0)  // user specified a number of events to generate
      params.generation().setMaxGen(num_events);

    if (params.generation().enabled() && !outputs.empty())
      for (const auto& output : outputs)
        gen.runParameters().addEventExporter(cepgen::EventExporterFactory::get().build(output));

    CG_LOG << gen.runParameters();  // user-friendly printout of all run parameters

    gen.computeXsection();  //--- let there be a cross-section...

    if (params.generation().enabled())
      // events generation happens here
      // (one may use a callback function)
      gen.generate(0);
  } catch (const cepgen::utils::RunAbortedException&) {
    CG_DEBUG("main") << "Run aborted!";
  } catch (const cepgen::Exception& e) {
    CG_DEBUG("main") << "CepGen exception encountered: " << e.what();
    e.dump();
  } catch (const exception& e) {
    CG_FATAL("main") << "Other exception caught!\n\t" << e.what();
  }
  return EXIT_SUCCESS;
}
