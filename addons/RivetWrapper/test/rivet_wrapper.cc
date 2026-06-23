#include <Rivet/Analysis.hh>
#include <Rivet/AnalysisHandler.hh>

#include "CepGen/EventFilter/EventExporter.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/EventExporterFactory.h"
#include "CepGen/Utils/Test.h"

using namespace std;

int main() {
  cepgen::initialise();

  vector<string> analyses = {"CMS_2011_I954992", "OPAL_1998_I474012"};

  auto* rivet_wrp = cepgen::EventExporterFactory::get()
                        .build("rivet", cepgen::ParametersList().set("analyses", analyses))
                        .release();  // do not call the destructor (Rivet will not be initialised)
  auto* rivet = rivet_wrp->engine<Rivet::AnalysisHandler>();
  CG_TEST_EQUAL(rivet->analysisNames(), analyses, "List of analyses");
  auto analysis = rivet->analysis(analyses.at(0));
  CG_TEST_EQUAL(analysis->experiment(), "CMS", "Analysis experiment");
  CG_TEST_EQUAL(analysis->collider(), "LHC", "Analysis collider");

  CG_TEST_SUMMARY;
}
