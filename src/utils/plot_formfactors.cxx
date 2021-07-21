#include "CepGen/Generator.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/FormFactors/Parameterisation.h"

#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

#include <fstream>

#include "TMultiGraph.h"
#include "TGraph.h"

using namespace std;

int main(int argc, char* argv[]) {
  int mode, strfun_type, num_points;
  double mx, q2min, q2max;
  string output_file;

  cepgen::ArgumentsParser(argc, argv)
      .addArgument("mode,t", "beam modelling", &mode, (int)cepgen::mode::Beam::ProtonElastic)
      .addOptionalArgument("mx,M", "diffractive mass (GeV/c^2)", &mx, 100.)
      .addOptionalArgument("sf,s", "structure functions modelling", &strfun_type, 301)
      .addOptionalArgument("q2min,m", "minimal parton virtuality (GeV^2)", &q2min, 1.)
      .addOptionalArgument("q2max,M", "maximal parton virtuality (GeV^2)", &q2max, 10000.)
      .addOptionalArgument("npoints,n", "number of x-points to scan", &num_points, 500)
      .addOptionalArgument("output,o", "output file name", &output_file, "formfacs.scan.output.txt")
      .parse();

  cepgen::initialise();
  //const double mi = cepgen::PDG::get().mass( cepgen::PDG::proton );
  //const double mi2 = mi*mi, mx2 = mx*mx;

  ofstream out(output_file);
  out << "# form factors: ";
  string sep;
  for (const auto& ff_type : cepgen::formfac::FormFactorsFactory::get().modules())
    out << sep << ff_type, sep = ", ";

  auto sf = cepgen::strfun::StructureFunctionsFactory::get().build(strfun_type);
  out << "\n"
      << "# structure functions: " << sf.get() << "\n"
      << "# q2 in [" << q2min << ", " << q2max << "] GeV^2\n";

  vector<unique_ptr<cepgen::formfac::Parameterisation> > form_factors;
  vector<TGraph*> g_form_factors_fe, g_form_factors_fm;
  for (const auto& ff_type : cepgen::formfac::FormFactorsFactory::get().modules()) {
    form_factors.emplace_back(cepgen::formfac::FormFactorsFactory::get().build(ff_type));
    (*form_factors.rbegin())->setStructureFunctions(sf.get());
    g_form_factors_fe.emplace_back(new TGraph);
    (*g_form_factors_fe.rbegin())->SetTitle((ff_type + ";Q^{2} (GeV^{2});F_{E}").c_str());
    g_form_factors_fm.emplace_back(new TGraph);
    (*g_form_factors_fm.rbegin())->SetTitle((ff_type + ";Q^{2} (GeV^{2});F_{M}").c_str());
  }
  for (int i = 0; i < num_points; ++i) {
    const double q2 = q2min + i * (q2max - q2min) / (num_points - 1);
    out << q2 << "\t";
    size_t j = 0;
    for (auto& ff : form_factors) {
      const auto form_factor = (*ff)((cepgen::mode::Beam)mode, q2, mx);
      out << "\t" << form_factor.FE << "\t" << form_factor.FM;
      g_form_factors_fe.at(j)->SetPoint(g_form_factors_fe.at(j)->GetN(), q2, form_factor.FE);
      g_form_factors_fm.at(j)->SetPoint(g_form_factors_fm.at(j)->GetN(), q2, form_factor.FM);
      ++j;
    }
    out << "\n";
  }
  CG_INFO("main") << "Scan written in \"" << output_file << "\".";
  out.close();

  for (auto& plt : map<const char*, vector<TGraph*> >{{"FE", g_form_factors_fe}, {"FM", g_form_factors_fm}}) {
    cepgen::ROOTCanvas c(plt.first, Form("M_{X} = %g GeV/c^{2}", mx));
    c.SetLogy();
    TMultiGraph mg;
    size_t i = 0;
    for (auto& gr : plt.second) {
      mg.Add(gr);
      gr->SetLineColor(cepgen::ROOTCanvas::colours[i]);
      c.AddLegendEntry(gr, gr->GetTitle(), "l");
      ++i;
    }
    mg.Draw("al");
    mg.GetHistogram()->SetTitle(Form(";Q^{2};%s", plt.first));
    c.Prettify(mg.GetHistogram());
    c.Save("pdf");
  }

  return 0;
}
