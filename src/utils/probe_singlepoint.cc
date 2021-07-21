#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main(int argc, char* argv[]) {
  string input_card;
  vector<double> point;
  bool enable_plugins, debug;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("input,i", "input card", &input_card)
      .addOptionalArgument("point,p", "point to test", &point, vector<double>(12, 0.3))
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .addOptionalArgument("enable-plugins,m", "enable the external plugins", &enable_plugins, false)
      .parse();

  cepgen::Generator gen;
  gen.setParameters(cepgen::card::Handler::parse(input_card));

  const auto ndim = gen.parameters()->process().ndim();
  if (point.size() < 2)
    point = vector<double>(ndim, point[0]);
  else if (point.size() != ndim)
    point.resize(ndim);

  if (debug)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debugInsideLoop;

  if (!enable_plugins) {
    gen.parametersPtr()->clearEventModifiersSequence();
    gen.parametersPtr()->clearOutputModulesSequence();
  }

  CG_INFO("main") << gen.parameters();

  CG_INFO("main") << "point: " << point;
  const double weight = gen.computePoint(point);
  CG_INFO("main") << "weight: " << weight;

  return 0;
}
