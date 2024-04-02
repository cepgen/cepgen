/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/AbortHandler.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Test.h"
#include "CepGen/Utils/Timer.h"

using namespace std;

int main(int argc, char* argv[]) {
  double num_sigma, chi2;
  int num_gen;
  string str_fun, proc_name, integrator, plotter;
  bool sublead_test;

  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("process,p", "process to compute", &proc_name, "lpair")
      .addOptionalArgument("num-gen,g", "number of events to generate", &num_gen, 50'000)
      .addOptionalArgument("num-sigma,n", "max. number of std.dev.", &num_sigma, 3.)
      .addOptionalArgument("str-fun,s", "struct.functions modelling", &str_fun, "SuriYennie")
      .addOptionalArgument("integrator,i", "type of integrator used", &integrator, "Vegas")
      .addOptionalArgument("plotter,t", "type of plotter to use", &plotter, "")
      .addOptionalArgument("chi2,x", "chi2 value cut for histograms compatibility test", &chi2, 1.)
      .addOptionalArgument("subleading-test", "also enable the subleading pt eta test?", &sublead_test, false)
      .parse();

  cepgen::utils::Timer tmr;
  cepgen::Generator gen;
  gen.runParameters().integrator().setName(integrator);

  cepgen::utils::AbortHandler ah;

  auto pkin =
      cepgen::ParametersList()
          .set<double>("sqrtS", 13.e3)
          .set<cepgen::ParametersList>(
              "structureFunctions", cepgen::StructureFunctionsFactory::get().describeParameters(str_fun).parameters())
          .set<double>("ptmin", 5.)
          .set<cepgen::Limits>("eta", {-2.5, 2.5})
          .set<cepgen::Limits>("mx", {1.07, 1000.});

  gen.runParameters().setProcess(
      cepgen::ProcessFactory::get().build(proc_name, cepgen::ParametersList().set<int>("pair", 13)));
  cepgen::Value cs_ei, cs_ie;

  auto h_eta_lead_ei = cepgen::utils::Hist1D(50, {-2.5, 2.5}, "eta_lead_ei", "el-inel"),
       h_eta_lead_ie = cepgen::utils::Hist1D(50, {-2.5, 2.5}, "eta_lead_ie", "inel-el"),
       h_eta_sublead_ei = cepgen::utils::Hist1D(50, {-2.5, 2.5}, "eta_sublead_ei", "el-inel"),
       h_eta_sublead_ie = cepgen::utils::Hist1D(50, {-2.5, 2.5}, "eta_sublead_ie", "inel-el"),
       h_mdiff_ei = cepgen::utils::Hist1D(50, {0., 1000.}, "mdiff_ei", "el-inel"),
       h_mdiff_ie = cepgen::utils::Hist1D(50, {0., 1000.}, "mdiff_ie", "inel-el");

  {  // elastic-inelastic
    pkin.set<int>("mode", 2);
    gen.runParameters().process().kinematics().setParameters(pkin);
    cs_ei = gen.computeXsection();
    if (num_gen > 0)
      gen.generate(num_gen, [&](const cepgen::Event& evt, size_t) {
        const auto &mom1 = evt(cepgen::Particle::Role::CentralSystem).at(0).momentum(),
                   &mom2 = evt(cepgen::Particle::Role::CentralSystem).at(1).momentum();
        if (mom1.pt() > mom2.pt()) {
          h_eta_lead_ei.fill(mom1.eta());
          h_eta_sublead_ei.fill(mom2.eta());
        } else {
          h_eta_lead_ei.fill(mom2.eta());
          h_eta_sublead_ei.fill(mom1.eta());
        }
        h_mdiff_ei.fill(evt(cepgen::Particle::Role::OutgoingBeam2).at(0).momentum().mass());
      });
  }
  {  // inelastic-elastic
    pkin.set<int>("mode", 3);
    gen.runParameters().process().kinematics().setParameters(pkin);
    cs_ie = gen.computeXsection();
    if (num_gen > 0)
      gen.generate(num_gen, [&](const cepgen::Event& evt, size_t) {
        const auto &mom1 = evt(cepgen::Particle::Role::CentralSystem).at(0).momentum(),
                   &mom2 = evt(cepgen::Particle::Role::CentralSystem).at(1).momentum();
        if (mom1.pt() > mom2.pt()) {
          h_eta_lead_ie.fill(mom1.eta());
          h_eta_sublead_ie.fill(mom2.eta());
        } else {
          h_eta_lead_ie.fill(mom2.eta());
          h_eta_sublead_ie.fill(mom1.eta());
        }
        h_mdiff_ie.fill(evt(cepgen::Particle::Role::OutgoingBeam1).at(0).momentum().mass());
      });
  }
  CG_TEST_VALUES(cs_ei, cs_ie, num_sigma, "el-inel == inel-el");

  size_t ndf;
  CG_TEST(h_eta_lead_ei.chi2test(h_eta_lead_ie, ndf) / ndf > chi2, "leading lepton eta");
  if (sublead_test)
    CG_TEST(h_eta_sublead_ei.chi2test(h_eta_sublead_ie, ndf) / ndf > chi2, "subleading lepton eta");
  CG_TEST(h_mdiff_ei.chi2test(h_mdiff_ie, ndf) / ndf < 1.5 * chi2, "diffractive system mass");

  if (!plotter.empty()) {
    auto plt = cepgen::DrawerFactory::get().build(plotter);
    plt->draw({&h_eta_lead_ie, &h_eta_lead_ei},
              "leading_eta",
              "leading lepton $\\eta$",
              cepgen::utils::Drawer::Mode::nostack);
    plt->draw({&h_eta_sublead_ie, &h_eta_sublead_ei},
              "subleading_eta",
              "subleading lepton $\\eta$",
              cepgen::utils::Drawer::Mode::nostack);
    plt->draw({&h_mdiff_ie, &h_mdiff_ei}, "mdiff", "diffractive system mass", cepgen::utils::Drawer::Mode::nostack);
  }

  CG_TEST_SUMMARY;
}
