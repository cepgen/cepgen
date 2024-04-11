/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

using namespace std::string_literals;

namespace cepgen {
  namespace apfelpp {
    class EvolutionStructureFunctions final : public strfun::Parameterisation {
    public:
      explicit EvolutionStructureFunctions(const ParametersList& params)
          : strfun::Parameterisation(params),
            apfel_grid_{{apfel::SubGrid{100, 1e-5, 3},
                         apfel::SubGrid{60, 1e-1, 3},
                         apfel::SubGrid{50, 6e-1, 3},
                         apfel::SubGrid{50, 8e-1, 3}}},
            proc_(steer<std::string>("processDIS")) {
        const auto masses = steer<std::vector<double> >("masses"),
                   thresholds = steer<std::vector<double> >("thresholds");
        const auto perturb_order = steer<int>("perturbativeOrder");
        const auto process_dis = steer<std::string>("processDIS");

        apfel::AlphaQCD alpha_s{0.35, sqrt(2), thresholds, perturb_order};
        const apfel::TabulateObject<double> alphas_tab{alpha_s, 100, 0.9, 1001, 3};
        const auto as = [&](double const& mu) -> double { return alphas_tab.Evaluate(mu); };

        // Effective charges
        auto fBq = [=](double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

        const auto dglap_obj = InitializeDglapObjectsQCD(apfel_grid_, thresholds);
        auto evolved_pdfs = apfel::BuildDglap(dglap_obj, apfel::LHToyPDFs, steer<double>("mu0"), perturb_order, as);
        const apfel::TabulateObject<apfel::Set<apfel::Distribution> > tabulated_pdfs{*evolved_pdfs, 50, 1, 1000, 3};
        const auto pdfs = [&](double const& x, double const& Q) -> std::map<int, double> {
          return tabulated_pdfs.EvaluateMapxQ(x, Q);
        };

        std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)> f2_obj, fl_obj;
        if (process_dis == "NC") {
          f2_obj = apfel::InitializeF2NCObjectsMassive(apfel_grid_, masses);
          fl_obj = apfel::InitializeFLNCObjectsMassive(apfel_grid_, masses);
        } else if (process_dis == "CC") {
          /*const auto F2PlusCCObj = InitializeF2CCPlusObjectsZM(apfel_grid_, thresholds);
          const auto FLPlusCCObj = InitializeFLCCPlusObjectsZM(apfel_grid_, thresholds);
          const auto F2MinusCCObj = InitializeF2CCMinusObjectsZM(apfel_grid_, thresholds);
          const auto FLMinusCCObj = InitializeFLCCMinusObjectsZM(apfel_grid_, thresholds);*/
        }
        f2_ = apfel::BuildStructureFunctions(f2_obj, pdfs, perturb_order, as, fBq);
        fl_ = apfel::BuildStructureFunctions(fl_obj, pdfs, perturb_order, as, fBq);

        // Tabulate Structure functions
        f2_total_.reset(new apfel::TabulateObject<apfel::Distribution>{
            [&](double const& Q) -> apfel::Distribution { return f2_.at(0).Evaluate(Q); }, 50, 1, 200, 3, thresholds});
        fl_total_.reset(new apfel::TabulateObject<apfel::Distribution>{
            [&](double const& Q) -> apfel::Distribution { return fl_.at(0).Evaluate(Q); }, 50, 1, 200, 3, thresholds});
      }

      static ParametersDescription description() {
        auto desc = strfun::Parameterisation::description();
        desc.setDescription("APFEL++ DIS structure functions");
        desc.add("mu0", M_SQRT2).setDescription("initial scale");
        desc.add("masses", std::vector<double>{0., 0., 0., M_SQRT2, 4.5, 175.});
        desc.add("thresholds", std::vector<double>{0., 0., 0.});
        desc.add("processDIS", "NC"s)
            .setDescription("process of the structure functions (NC, or CC)")
            .allow("NC", "neutral currents")
            .allow("CC", "charged currents");
        return desc;
      }

    private:
      void eval() override {
        setF2(f2_total_->EvaluatexQ(args_.xbj, args_.q2));
        setFL(fl_total_->EvaluatexQ(args_.xbj, args_.q2));
      }

      const apfel::Grid apfel_grid_;  // x-space grid
      const std::string proc_;

      std::map<int, apfel::DglapObjects> dglap_obj_;

      std::map<int, apfel::Observable<> > f2_, fl_;
      std::unique_ptr<apfel::TabulateObject<apfel::Distribution> > f2_total_, fl_total_;
    };
  }  // namespace apfelpp
}  // namespace cepgen
using EvolutionStructureFunctionsAPFELpp = cepgen::apfelpp::EvolutionStructureFunctions;
REGISTER_STRFUN("apfelppEvol", 405, EvolutionStructureFunctionsAPFELpp);
