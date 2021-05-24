#include "CepGen/Utils/ArgumentsParser.h"

#include <LHAPDF/LHAPDF.h>

using namespace std;

int main(int argc, char* argv[]) {
  double q2;
  string set;
  int member;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("q2", "Virtuality", &q2, 100.)
      .addOptionalArgument("set,s", "PDFset to use", &set, "LUXqed17_plus_PDF4LHC15_nnlo_100")
      .addOptionalArgument("member,m", "PDF member", &member, 0)
      .parse();

  const LHAPDF::PDF* pdf = LHAPDF::mkPDF(set, member);
  double x = 0.;
  const double xfx = pdf->xfxQ2(22, x, q2);
  return 0;
}
