#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/ArgumentsParser.h"

using namespace std;

int main(int argc, char* argv[]) {
  string path;
  cepgen::ArgumentsParser parser(argc, argv);
  parser.addOptionalArgument("input,i", "path to the MCD file", &path, "../External/mass_width_2019.mcd").parse();
  pdg::MCDFileParser::parse(path);
  cepgen::PDG::get().dump();
  return 0;
}
