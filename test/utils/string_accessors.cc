#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"

using namespace std;

int main(int argc, char* argv[]) {
  cepgen::ArgumentsParser(argc, argv).parse();

  {
    auto str = "This {is} a text with {sub-strings}";
    const auto substr = cepgen::utils::between(str, "{", "}");
    if (substr.size() != 2) {
      CG_LOG << "Invalid number of substring found: " << substr << ".";
      return -1;
    }
  }

  return 0;
}
