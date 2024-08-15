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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/EventFilter/EventImporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/EventImporterFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_file, output_file;

  cepgen::ArgumentsParser parser(argc, argv);
  parser.addArgument("input,i", "input event file", &input_file)
      .addArgument("output,o", "output event file", &output_file)
      .parse();

  cepgen::initialise();

  auto params = cepgen::RunParameters{};

  auto reader = cepgen::EventImporterFactory::get().build(input_file);
  reader->initialise(params);
  auto writer = cepgen::EventExporterFactory::get().build(output_file);
  writer->initialise(params);

  writer->setCrossSection(reader->crossSection());

  cepgen::Event buf;
  size_t num_events_converted = 0;
  while ((*reader) >> buf) {
    (*writer) << buf;
    ++num_events_converted;
  }

  CG_LOG << "Successfully converted " << cepgen::utils::s("event", num_events_converted, true) << ".";

  return 0;
}
