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

/// Photon-photon (UPC-like, elastic-elastic) validation for a configurable set of initial states.
///
/// For each initial state (pp, pPb, PbPb, ep, see InitialStates.h), and for each final state (lepton pair or W pair),
/// the distributions of the central system are compared between:
///  * LPAIR (lepton pairs only),
///  * the kt-factorised processes (PPtoFF, PPtoWW), for one or several families of kt-dependent photon flux modellings.
/// One set of plots is produced per (final state, initial state) combination.

#include "CepGen/Core/RunParameters.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/DrawerFactory.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/Process.h"
#include "CepGen/Utils/ArgumentsParser.h"
#include "CepGen/Utils/Collections.h"
#include "CepGen/Utils/Drawer.h"
#include "CepGen/Utils/Histogram.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "Comparator.h"
#include "InitialStates.h"

using namespace std;
using namespace std::string_literals;

namespace {
  /// Pair-level histograms for a given initial/final state combination
  struct comparison_t : public cepgen::validation::Comparator {
    comparison_t(cepgen::Generator& gen,
                 const cepgen::ParametersList& params,
                 const cepgen::validation::InitialState& state,
                 bool ww)
        : cepgen::validation::Comparator(gen, params), state_(state), ww_(ww) {}

    void initialise() override {
      const auto mass_range = ww_ ? cepgen::Limits{150., state_.ww_mass_max} : cepgen::Limits{0., state_.mass_max};
      const auto pt_max = ww_ ? 50. : state_.pt_pair_max;
      const auto pair = ww_ ? "W^{+}W^{-}"s : "l^{+}l^{-}"s;
      (*this)
          .book("invmass", "$m(" + pair + ")$", "GeV", cepgen::utils::Hist1D(50, mass_range))
          .book("ptpair", "$p_{T}(" + pair + ")$", "GeV", cepgen::utils::Hist1D(50, {0., pt_max}))
          .book("ypair", "$y(" + pair + ")$", "", cepgen::utils::Hist1D(50, {-6., 6.}))
          .book("ptlead", "$p_{T}^{lead}$", "GeV", cepgen::utils::Hist1D(50, {0., ww_ ? 250. : 4. * state_.ptmin}))
          .book("etalead", "$\\eta^{lead}$", "", cepgen::utils::Hist1D(50, {-6., 6.}))
          .book("etasublead", "$\\eta^{sublead}$", "", cepgen::utils::Hist1D(50, {-6., 6.}))
          .book("acop", "$1-|\\Delta\\phi/\\pi|$", "", cepgen::utils::Hist1D(50, {0., ww_ ? 1. : 0.5}));
      for (const auto& plot : {"invmass", "ptpair", "ptlead", "acop"})
        drawMode(plot) |= cepgen::utils::Drawer::Mode::logy;
    }
    void process(const cepgen::Event& evt) override {
      const auto &cm = evt(cepgen::Particle::Role::Intermediate).at(0).momentum(),
                 &pl1 = evt(cepgen::Particle::Role::CentralSystem).at(0).momentum(),
                 &pl2 = evt(cepgen::Particle::Role::CentralSystem).at(1).momentum();
      const auto &lead = pl1.pt() > pl2.pt() ? pl1 : pl2, &sublead = pl1.pt() > pl2.pt() ? pl2 : pl1;
      (*this)
          .fill("invmass", cm.mass())
          .fill("ptpair", cm.pt())
          .fill("ypair", cm.rapidity())
          .fill("ptlead", lead.pt())
          .fill("etalead", lead.eta())
          .fill("etasublead", sublead.eta())
          .fill("acop", 1. - fabs(pl1.deltaPhi(pl2) * M_1_PI));
    }

  private:
    const cepgen::validation::InitialState& state_;
    const bool ww_;
  };
}  // namespace

int main(int argc, char* argv[]) {
  int num_gen, num_threads;
  double eta_max;
  vector<string> states, processes, flux_families;
  string filename, plotter, final_state;
  cepgen::ArgumentsParser(argc, argv)
      .addOptionalArgument("states,s",
                           "initial states to validate (any of: pp, pPb, PbPb, ep)",
                           &states,
                           vector{"pp"s, "pPb"s, "PbPb"s, "ep"s})
      .addOptionalArgument("final-state,f", "final state: 'll' (lepton pair) or 'ww'", &final_state, "ll"s)
      .addOptionalArgument("processes,P",
                           "processes to generate (lpair is only used for 'll', pptoww for 'ww')",
                           &processes,
                           vector{"lpair"s, "pptoff"s, "pptoww"s})
      .addOptionalArgument("fluxes,F",
                           "kt-factorised flux modelling families (any of: default, elastic, budnev, klein)",
                           &flux_families,
                           vector{"default"s})
      .addOptionalArgument("eta-max,e",
                           "single-lepton |eta| cut (<= 0 to disable; useful to check frame-dependent effects)",
                           &eta_max,
                           2.5)
      .addOptionalArgument("num-gen,n", "number of events to generate", &num_gen, 10'000)
      .addOptionalArgument("num-threads,t", "number of threads to generate events", &num_threads, 4)
      .addOptionalArgument("plotter,p", "type of plotter to user", &plotter, "root")
      .addOptionalArgument(
          "filename,o", "output base filename", &filename, "validation/comparison_pair_initial_states_")
      .parse();

  const bool ww = (final_state == "ww");
  if (!ww && final_state != "ll")
    throw CG_FATAL("main") << "Invalid final state: '" << final_state << "'. Allowed values: ll, ww.";

  cepgen::Generator gen;
  for (const auto& state_name : states) {
    const auto& state = cepgen::validation::InitialState::get(state_name);
    comparison_t comp(
        gen,
        cepgen::ParametersList()
            .set("topLabel",
                 state.label + ", elastic-elastic $\\gamma\\gamma \\rightarrow "s + (ww ? "W^{+}W^{-}" : "l^{+}l^{-}") +
                     "$")
            .set("numEvents", num_gen)
            .set("numThreads", num_threads)
            .set("pathTemplate", filename + final_state + "_" + state.name + (eta_max > 0. ? "" : "_noetacut") + "_")
            .set("plotter", cepgen::ParametersList().setName(plotter).feed(plotter).set("format", "png,pdf"s)),
        state,
        ww);

    auto& pars = gen.runParameters();
    for (const auto& proc_name : processes) {
      if (!cepgen::utils::contains(cepgen::ProcessFactory::get().modules(), proc_name)) {
        CG_WARNING("main") << "Process '" << proc_name << "' is not available. Skipping.";
        continue;
      }
      const bool is_lpair = (proc_name == "lpair");
      if ((ww && proc_name != "pptoww") || (!ww && proc_name == "pptoww"))
        continue;  // process incompatible with the requested final state

      // LPAIR has a built-in photon flux: only one sample is generated for it
      for (const auto& flux : is_lpair ? vector{"default"s} : flux_families) {
        auto proc_params = cepgen::ParametersList();
        if (!is_lpair)
          proc_params.set<std::string>("partonsGenerator", "kt");
        if (ww)
          proc_params.set<int>("method", 1 /* off-shell */).set<int>("polarisationStates", 0 /* full */);
        else
          proc_params.set<int>("pair", 13 /* muons */);
        pars.setProcess(cepgen::ProcessFactory::get().build(proc_name, proc_params));

        auto kin_params = state.kinematics(is_lpair ? "default"s : flux);
        if (!ww) {
          kin_params.set<double>("ptmin", state.ptmin);
          if (eta_max > 0.)
            kin_params.set<cepgen::Limits>("eta", {-eta_max, eta_max});
        }
        pars.process().kinematics().setParameters(kin_params);
        comp.loop(proc_name + (is_lpair || flux == "default" ? ""s : " (" + flux + ")"));
      }
    }
  }
  return 0;
}
