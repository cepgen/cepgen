/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2025  Laurent Forthomme
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

#include <apfel/apfelxx.h>

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/CouplingFactory.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/Coupling.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

using namespace std::string_literals;

namespace cepgen::apfelpp {
  class EvolutionStructureFunctions final : public strfun::Parameterisation {
  public:
    explicit EvolutionStructureFunctions(const ParametersList& params)
        : Parameterisation(params),
          apfel_grid_{{apfel::SubGrid{100, 1e-5, 3},
                       apfel::SubGrid{60, 1e-1, 3},
                       apfel::SubGrid{50, 6e-1, 3},
                       apfel::SubGrid{50, 8e-1, 3}}},
          alpha_s_(AlphaSFactory::get().build(steer<ParametersList>("alphaSParameters"))) {
      if (alpha_s_->name() != "apfelpp")  // do not print twice the APFEL++ banner
        apfel::Banner();
      const auto thresholds = steer<std::vector<double> >("thresholds");
      const auto mu0 = steer<double>("mu0");
      const auto perturb_order = steer<int>("perturbativeOrder");

      const auto as = [&](double mu) { return (*alpha_s_)(mu); };
      const auto fBq = [=](double q) { return apfel::ElectroWeakCharges(q, false); };  // effective charges

      const auto dglap_obj = InitializeDglapObjectsQCD(apfel_grid_, thresholds);
      auto evolved_pdfs = BuildDglap(dglap_obj, apfel::LHToyPDFs, mu0, perturb_order, as);
      const apfel::TabulateObject tabulated_pdfs{*evolved_pdfs, 50, 1, 1000, 3};
      const auto pdfs = [&](double x, double q) { return tabulated_pdfs.EvaluateMapxQ(x, q); };

      const auto process_dis = steer<std::string>("processDIS");
      if (process_dis == "NC") {
        const auto masses = steer<std::vector<double> >("masses");
        const auto f2_obj = InitializeF2NCObjectsMassive(apfel_grid_, masses);
        const auto fl_obj = InitializeFLNCObjectsMassive(apfel_grid_, masses);
        auto f2 = BuildStructureFunctions(f2_obj, pdfs, perturb_order, as, fBq);
        auto fl = BuildStructureFunctions(fl_obj, pdfs, perturb_order, as, fBq);
        // tabulate structure functions
        f2_total_.reset(new apfel::TabulateObject<apfel::Distribution>{
            [&](double q) { return f2.at(0).Evaluate(q); }, 50, 1, 200, 3, thresholds});
        fl_total_.reset(new apfel::TabulateObject<apfel::Distribution>{
            [&](double q) { return fl.at(0).Evaluate(q); }, 50, 1, 200, 3, thresholds});
      } else if (process_dis == "CC") {
        const auto f2p_obj = InitializeF2CCPlusObjectsZM(apfel_grid_, thresholds);
        const auto f2m_obj = InitializeF2CCMinusObjectsZM(apfel_grid_, thresholds);
        const auto flp_obj = InitializeFLCCPlusObjectsZM(apfel_grid_, thresholds);
        const auto flm_obj = InitializeFLCCMinusObjectsZM(apfel_grid_, thresholds);
        auto f2p = BuildStructureFunctions(f2p_obj, pdfs, perturb_order, as, fBq);
        auto f2m = BuildStructureFunctions(f2m_obj, pdfs, perturb_order, as, fBq);
        auto flp = BuildStructureFunctions(flp_obj, pdfs, perturb_order, as, fBq);
        auto flm = BuildStructureFunctions(flm_obj, pdfs, perturb_order, as, fBq);
        // tabulate structure functions
        f2_total_.reset(new apfel::TabulateObject<apfel::Distribution>{
            [&](double q) { return f2p.at(0).Evaluate(q) - f2m.at(0).Evaluate(q); }, 50, 1, 200, 3, thresholds});
        fl_total_.reset(new apfel::TabulateObject<apfel::Distribution>{
            [&](double q) { return flp.at(0).Evaluate(q) - flm.at(0).Evaluate(q); }, 50, 1, 200, 3, thresholds});
      }
    }

    static ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("APFEL++ DIS structure functions");
      desc.add("mu0", 1.3).setDescription("initial scale");
      desc.add("masses", std::vector<double>{0., 0., 0., 1.3, 4.75, 175.});
      desc.add("thresholds", std::vector<double>{0., 0., 0.});
      desc.add("perturbativeOrder", 0)
          .setDescription("perturbative order for alpha(S) evolution")
          .allow(0, "LO")
          .allow(1, "NLO")
          .allow(2, "NNLO")
          .allow(3, "NNNLO");
      desc.add("processDIS", "NC"s)
          .setDescription("process of the structure functions (NC, or CC)")
          .allow("NC", "neutral currents")
          .allow("CC", "charged currents");
      desc.add("alphaSParameters", AlphaSFactory::get().describeParameters("apfelpp"));
      return desc;
    }

  private:
    void eval() override {
      const auto q = std::sqrt(args_.q2);
      setF2(f2_total_->EvaluatexQ(args_.xbj, q));
      setFL(fl_total_->EvaluatexQ(args_.xbj, q));
    }

    const apfel::Grid apfel_grid_;  // x-space grid
    const std::unique_ptr<Coupling> alpha_s_;
    std::unique_ptr<apfel::TabulateObject<apfel::Distribution> > f2_total_, fl_total_;
  };
}  // namespace cepgen::apfelpp
using EvolutionStructureFunctionsAPFELpp = cepgen::apfelpp::EvolutionStructureFunctions;
REGISTER_STRFUN("apfelppEvol", 405, EvolutionStructureFunctionsAPFELpp);
