#include <Rivet/Analysis.hh>
#include <Rivet/AnalysisHandler.hh>

#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  cepgen::Generator gen;

  auto proc = cepgen::ProcessFactory::get().build("lpair");
  gen.runParameters().setProcess(std::move(proc));

  auto& kin = gen.runParameters().process().kinematics();
  kin.incomingBeams().positive().setIntegerPdgId(2212);
  kin.incomingBeams().negative().setIntegerPdgId(2212);
  kin.incomingBeams().setSqrtS(7.e3);

  auto rivet_wrp = cepgen::EventExporterFactory::get().build(
      "rivet", cepgen::ParametersList().set<std::vector<std::string> >("analyses", {"CMS_2011_I954992"}));
  //rivet_wrp->initialise(gen.runParameters());
  gen.runParameters().addEventExporter(std::move(rivet_wrp));

  gen.generate(100);

  return 0;
}
