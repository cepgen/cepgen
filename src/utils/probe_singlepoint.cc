#include "CepGen/Cards/Handler.h"
#include "CepGen/Generator.h"
#include "CepGen/Parameters.h"
#include "CepGen/Processes/Process.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main(int argc, char* argv[]) {
  const size_t ps_size = 12;
  string input_card;
  vector<double> point;
  bool enable_plugins, debug;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("input,i", "input card", &input_card)
      .addOptionalArgument("point,p", "point to test", &point, vector<double>(ps_size, 0.3))
      .addOptionalArgument("debug,d", "debugging mode", &debug, false)
      .addOptionalArgument("enable-plugins,m", "enable the external plugins", &enable_plugins, false)
      .parse();

  cepgen::Generator gen;
  gen.setParameters(cepgen::card::Handler::parse(input_card));

  const auto ndim = gen.parameters()->process().ndim();
  if (point.size() < 2) {
    point = vector<double>(ndim, point[0]);
    point.resize(ps_size);
  } else if (point.size() != ndim)
    point.resize(ndim);

  if (!enable_plugins) {
    gen.parametersPtr()->clearEventModifiersSequence();
    gen.parametersPtr()->clearOutputModulesSequence();
  }

  CG_INFO("main") << gen.parameters();

  if (debug)
    cepgen::utils::Logger::get().level = cepgen::utils::Logger::Level::debugInsideLoop;

  cout << "point: ";
  string delim;
  for (const auto& v : point)
    cout << delim << v, delim = ", ";
  cout << endl;
  const double weight = gen.computePoint(point);
  cout << "weight: " << weight << endl;

  return 0;
}
