/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#include <cmath>

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Physics/Modes.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "Comparator.h"

using namespace std;
using namespace std::string_literals;

/// Human-readable label for a CepGen kinematics mode (cepgen::mode::Kinematics numbering)
static std::string modeFromId(int mode) {
  switch (static_cast<cepgen::mode::Kinematics>(mode)) {
    case cepgen::mode::Kinematics::ElasticElastic:
      return "Elastic";
    case cepgen::mode::Kinematics::InelasticElastic:  // beam 1 dissociates
      return "SD^{+}";
    case cepgen::mode::Kinematics::ElasticInelastic:  // beam 2 dissociates
      return "SD^{-}";
    case cepgen::mode::Kinematics::InelasticInelastic:
      return "DD";
    default:
      return "invalid";
  }
}

/// Compare the kT-factorised two-photon production of a pair of charged scalars (sleptons) between CepGen's pptoss and
/// mg5_aMC, using the MSSM_SLHA2 model for the latter.
/// \note The slepton is defined with the same PDG id and mass in both generators: the mass of the MadGraph particle is
///  steered through the MASS block of the model parameters card.
/// \note In MSSM_SLHA2, the slepton mass eigenstates follow the SLHA2 numbering (PDG ids 1000011, 1000013, 1000015,
///  2000011, 2000013, 2000015), and the MadGraph particle name must be consistent with the chosen PDG id. Run
///  "import model MSSM_SLHA2" then "display particles" in mg5_aMC to check the name matching --pdgid, and set it with
///  --mg5-process otherwise.
int main(int argc, char* argv[]) {
  int num_gen, num_threads, mode, pdg_id;
  double mass, sqrt_s, pt_min, qt_max;
  vector<string> processes;
  string filename, plotter, mg5_process, mg5_model;
  bool ratio_plot;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("processes,P", "processes to generate", &processes, vector<string>{"pptoss", "mg5_aMC"})
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 10'000)
      .addOptionalArgument("num-threads,t", "number of threads to generate events", &num_threads, 4)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "root")
      .addOptionalArgument("ratio,r", "draw the ratio plot", &ratio_plot, false)
      .addOptionalArgument("filename,f", "output base filename", &filename, "validation/comparison_slepton_kt_")
      .addOptionalArgument("mode", "kinematics mode (1=el-el, 2=el-inel, 3=inel-el, 4=inel-inel)", &mode, 3)
      .addOptionalArgument("pdgid", "slepton PDG id (SLHA2 numbering)", &pdg_id, 1000011)
      .addOptionalArgument("mass,m", "slepton mass (in GeV)", &mass, 100.)
      .addOptionalArgument("sqrts,s", "centre of mass energy (in GeV)", &sqrt_s, 13.6e3)
      .addOptionalArgument("ptmin", "minimum slepton pT (in GeV)", &pt_min, 10.)
      .addOptionalArgument("qtmax", "maximum parton kT (in GeV)", &qt_max, 50.)
      .addOptionalArgument("mg5-model", "mg5_aMC model", &mg5_model, "MSSM_SLHA2")
      .addOptionalArgument("mg5-process", "mg5_aMC process", &mg5_process, "a a > el- el+")
      .parse();

  // define the slepton for CepGen's own PDG library (charge in units of e/3)
  cepgen::PDG::get().define(
      cepgen::ParticleProperties(static_cast<cepgen::pdgid_t>(pdg_id), "slepton", "slepton", 0, mass, 0., {-3}, false));

  struct comparison_t : public cepgen::validation::Comparator {
    using cepgen::validation::Comparator::Comparator;
    void initialise() override {
      (*this)
          .book("invmass", "$m(\\tilde{l}^{+}\\tilde{l}^{-})$", "GeV", cepgen::utils::Hist1D(72, {200., 2000.}))
          .book("ptpair", "$p_{T}(\\tilde{l}^{+}\\tilde{l}^{-})$", "GeV", cepgen::utils::Hist1D(50, {0., 50.}))
          .book("ptlead", "$p_{T}^{lead}$", "GeV", cepgen::utils::Hist1D(50, {0., 200.}))
          .book("ptsublead", "$p_{T}^{sublead}$", "GeV", cepgen::utils::Hist1D(50, {0., 200.}))
          .book("etalead", "$\\eta^{lead}$", "", cepgen::utils::Hist1D(50, {-5., 5.}))
          .book("etasublead", "$\\eta^{sublead}$", "", cepgen::utils::Hist1D(50, {-5., 5.}))
          .book(
              "acop", "$1-|\\Delta\\phi(\\tilde{l}^{+}\\tilde{l}^{-})/\\pi|$", "", cepgen::utils::Hist1D(50, {0., 1.}))
          .book("mx", "$M_{X}$", "GeV", cepgen::utils::Hist1D(50, {0., 1000.}));
      for (const auto& plot : {"invmass", "ptpair", "ptlead", "ptsublead", "acop", "mx"})
        drawMode(plot) |= cepgen::utils::Drawer::Mode::logy;
    }
    void process(const cepgen::Event& evt) override {
      const auto &cm = evt(cepgen::Particle::Role::Intermediate).at(0).momentum(),
                 &pl1 = evt(cepgen::Particle::Role::CentralSystem).at(0).momentum(),
                 &pl2 = evt(cepgen::Particle::Role::CentralSystem).at(1).momentum();
      cepgen::Momentum pl_lead, pl_sublead;
      if (pl1.pt() > pl2.pt())
        pl_lead = pl1, pl_sublead = pl2;
      else
        pl_lead = pl2, pl_sublead = pl1;
      (*this)
          .fill("invmass", cm.mass())
          .fill("ptpair", cm.pt())
          .fill("ptlead", pl_lead.pt())
          .fill("etalead", pl_lead.eta())
          .fill("ptsublead", pl_sublead.pt())
          .fill("etasublead", pl_sublead.eta())
          .fill("acop", 1. - fabs(pl1.deltaPhi(pl2) * M_1_PI))
          .fill("mx", evt(cepgen::Particle::Role::OutgoingBeam1).at(0).momentum().mass());
    }
  };

  cepgen::Generator gen;
  auto comp = comparison_t(
      gen,
      cepgen::ParametersList()
          .set(
              "topLabel",
              modeFromId(mode) +
                  cepgen::utils::format(
                      " $\\gamma\\gamma \\rightarrow \\tilde{l}^{+}\\tilde{l}^{-}$ (%g TeV), $p_{T}^{\\tilde{l}} >$ %g GeV, $k_{T}$"s,
                      sqrt_s * 1.e-3,
                      pt_min))
          .set("numEvents", num_gen)
          .set("numThreads", num_threads)
          .set("pathTemplate", filename)
          .set("plotter", cepgen::ParametersList().setName(plotter).feed(plotter).set("format", "png,pdf"s)));

  auto& pars = gen.runParameters();
  for (const auto& proc_name : processes) {
    const auto module_name = cepgen::utils::split(proc_name, '<').at(0);
    auto proc = proc_name;
    if (!cepgen::utils::contains(cepgen::ProcessFactory::get().modules(), module_name)) {
      CG_WARNING("compare_slepton_kt") << "Process module '" << module_name << "' is not registered; skipping.";
      continue;
    }
    auto proc_params = cepgen::ParametersList().set("partonsGenerator", "kt"s);
    if (module_name == "mg5_aMC") {
      proc += "<process:'"s + mg5_process + "'";
      proc += "<model:'"s + mg5_model + "'";
      proc += "<modelParameters/MASS/"s + std::to_string(pdg_id) + ":" + std::to_string(mass);  // slepton mass (SLHA)
    } else
      proc_params.set<int>("pair", pdg_id);
    pars.setProcess(cepgen::ProcessFactory::get().build(proc, proc_params));
    pars.process().kinematics().setParameters(cepgen::ParametersList()
                                                  .set<vector<int> >("pdgIds", {2212, 2212})
                                                  .set<double>("sqrtS", sqrt_s)
                                                  .set<int>("mode", mode)
                                                  .set<double>("ptmin", pt_min)
                                                  .set<double>("qtmax", qt_max));
    comp.loop(proc_name);
  }
  return 0;
}
