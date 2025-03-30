#include "CepGen/Event/Event.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/RandomGeneratorFactory.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Filesystem.h"
#include "CepGen/Utils/RandomGenerator.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main(int argc, char* argv[]) {
  string output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("output,o", "path to the output LHEF file", &output_file, "test.lhe")
      .parse();

  cepgen::initialise();

  // initialise the LHEF writer
  auto lhef_mod =
      cepgen::EventExporterFactory::get().build("lhef", cepgen::ParametersList().set("filename", output_file));

  // randomise the number of events to be written in the output file
  const auto rng = cepgen::RandomGeneratorFactory::get().build("stl");
  const size_t num_events = rng->uniformInt(1, 10);

  // generate one simple event
  const auto evt = cepgen::Event::minimal();
  {  // write a few events
    for (size_t i = 0; i < num_events; ++i)
      (*lhef_mod) << evt;
  }

  // start of tests on output file
  CG_TEST_EQUAL(cepgen::utils::fileExists(output_file), true, "Output file exists");

  size_t num_stored_events = 0, num_events_invalid_multiplicity_hdr = 0, num_events_invalid_multiplicity_cnt = 0;

  bool in_event = false;
  size_t num_lines_in_event = 999, num_particles_in_event_hdr = 999;
  for (const auto& buf : cepgen::utils::split(cepgen::utils::readFile(output_file), '\n')) {
    if (buf == "<event>") {
      in_event = true;
      num_lines_in_event = 0;
    } else if (in_event && buf == "</event>") {
      if (num_lines_in_event - 2 != num_particles_in_event_hdr)  // skip incoming beam particles
        ++num_events_invalid_multiplicity_cnt;
      num_stored_events++;
      in_event = false;
    }
    if (num_lines_in_event == 1) {
      num_particles_in_event_hdr = stoul(cepgen::utils::split(cepgen::utils::trim(buf), ' ', true).at(0));
      if (num_particles_in_event_hdr - 2 == evt.size())  // remove </event> and header line
        ++num_events_invalid_multiplicity_hdr;
    }
    ++num_lines_in_event;
  }

  CG_TEST_EQUAL(num_stored_events, num_events, "Number of <event> + </event> tags in LHEF");
  CG_TEST_EQUAL(num_events_invalid_multiplicity_hdr,
                0,
                "No events with invalid header-registered particles multiplicity in LHEF");
  CG_TEST_EQUAL(num_events_invalid_multiplicity_cnt, 0, "No events with invalid particles multiplicity in LHEF");

  CG_TEST_SUMMARY;
}
