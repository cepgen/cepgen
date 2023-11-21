#include <Rivet/Analysis.hh>
#include <Rivet/AnalysisHandler.hh>

#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Parameters.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::Generator gen;

  auto proc = cepgen::ProcessFactory::get().build("lpair");
  gen.parametersRef().setProcess(std::move(proc));

  auto& kin = gen.parametersRef().process().kinematics();
  kin.incomingBeams().positive().setPdgId(2212);
  kin.incomingBeams().negative().setPdgId(2212);
  kin.incomingBeams().setSqrtS(7.e3);

  auto rivet_wrp = cepgen::EventExporterFactory::get().build(
      "rivet", cepgen::ParametersList().set<std::vector<std::string> >("analyses", {"CMS_2011_I954992"}));
  //rivet_wrp->initialise(gen.parametersRef());
  gen.parametersRef().addEventExporter(std::move(rivet_wrp));

  gen.generate();

  return 0;
}
