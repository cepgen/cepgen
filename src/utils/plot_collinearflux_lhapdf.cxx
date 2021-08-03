/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#include <LHAPDF/LHAPDF.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/CollinearFlux.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace std;

int main(int argc, char* argv[]) {
  double q2, xmin, xmax;
  string ffmode, set, output;
  int strfun_type, member, num_points;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("q2", "Virtuality", &q2, 100.)
      .addOptionalArgument("xmin,x", "minimal fractional loss", &xmin, 1.e-5)
      .addOptionalArgument("xmax,X", "maximal fractional loss", &xmax, 1.)
      .addOptionalArgument("formfac,f", "form factors modelling", &ffmode, "StandardDipole")
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("set,s", "PDFset to use", &set, "LUXqed17_plus_PDF4LHC15_nnlo_100")
      .addOptionalArgument("output,o", "Output filename", &output, argv[0])
      .addOptionalArgument("member,m", "PDF member", &member, 0)
      .addOptionalArgument("num-points,n", "Number of points to probe", &num_points, 100)
      .parse();

  cepgen::initialise();

  unique_ptr<LHAPDF::PDF> pdf(LHAPDF::mkPDF(set, member));

  //auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(
      401, cepgen::ParametersList().set<std::string>("pdfSet", set).set<int>("pdfMember", member));
  auto ff = cepgen::formfac::FormFactorsFactory::get().build(ffmode);
  ff->setStructureFunctions(sf.get());

  const cepgen::Limits kt2_limits(0., 1000.);

  const cepgen::CollinearFlux flux(ff.get(), kt2_limits);

  TGraph g_ref, g_cg, g_ratio;
  for (int i = 0; i < num_points; ++i) {
    const double x = xmin + i * (xmax - xmin) / (num_points + 1);
    const double xfx = pdf->xfxQ2(22, x, q2);
    const double pdf = flux(x, 0.938, cepgen::KTFlux::P_Photon_Elastic_Budnev);
    cout << x << "\t" << xfx << "\t" << pdf << "\t" << pdf / xfx << endl;
    g_ref.SetPoint(g_ref.GetN(), x, xfx);
    g_cg.SetPoint(g_cg.GetN(), x, pdf);
    g_ratio.SetPoint(g_ratio.GetN(), x, pdf / xfx);
  }

  cepgen::ROOTCanvas c(output.c_str());
  TMultiGraph mg;
  g_ref.SetLineColor(kRed + 1);
  g_cg.SetLineColor(kBlue + 2);
  //mg.Add(&g_ref);
  //mg.Add(&g_cg);
  mg.Add(&g_ratio);
  mg.SetMinimum(1.e-10);
  mg.Draw("al");
  c.Prettify(mg.GetHistogram());
  c.SetLogy();
  c.Save("pdf");
  return 0;
}
