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

#include <APFEL/APFEL.h>

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

using namespace std::string_literals;

namespace cepgen::apfel {
  class EvolutionStructureFunctions final : public strfun::Parameterisation {
  public:
    explicit EvolutionStructureFunctions(const ParametersList& params)
        : Parameterisation(params),
          proc_(steer<std::string>("processDIS")),
          q2_range_(steer<Limits>("q2range")),
          xbj_min_(steer<double>("xBjmin")) {
      if (!q2_range_.valid())
        throw CG_FATAL("apfel:EvolutionStructureFunctions") << "Invalid Q^2 range: " << q2_range_ << ".";
      const auto qrange = q2_range_.compute([](double lim) { return std::sqrt(lim); });
      APFEL::SetMassScheme(steer<std::string>("massScheme"));
      APFEL::SetProcessDIS(proc_);
      APFEL::SetQLimits(qrange.min(), qrange.max());
      APFEL::SetPerturbativeOrder(steer<int>("perturbativeOrder"));
      APFEL::SetMaxFlavourAlpha(steer<int>("maxFlavourAlpha"));
      APFEL::SetMaxFlavourPDFs(steer<int>("maxFlavourPDFs"));
      APFEL::SetPDFSet(steer<std::string>("pdfSet"));
      APFEL::SetTargetDIS(steer<std::string>("targetDIS"));
      APFEL::InitializeAPFEL_DIS();
      APFEL::ComputeStructureFunctionsAPFEL(qrange.min(), qrange.max());
      APFEL::CacheStructureFunctionsAPFEL(qrange.min());
    }

    static ParametersDescription description() {
      auto desc = strfun::Parameterisation::description();
      desc.setDescription("APFEL DIS structure functions");
      desc.add("q2range", Limits{1., 1.e6}).setDescription("evolution scale range, in GeV^2");
      desc.add<double>("xBjmin", 2.e-6).setDescription("minimum Bjorken-x reachable for this PDF set");
      desc.add("massScheme", "FFNS"s /*"FONLL-A"s*/).setDescription("mass scheme for the structure functions");
      desc.add("processDIS", "NC"s).setDescription("process of the structure functions (EM, NC, or CC)");
      desc.add("perturbativeOrder", 0)
          .setDescription("perturbative order for alpha(S) evolution")
          .allow(0, "LO")
          .allow(1, "NLO")
          .allow(2, "NNLO")
          .allow(3, "NNNLO");
      desc.add<int>("maxFlavourAlpha", 5)
          .setDescription("maximum number of flavours that the evolution of alpha(S) and alpha(EM) can reach");
      desc.add<int>("maxFlavourPDFs", 5)
          .setDescription("maximum number of flavours that the evolution of PDFs can reach");
      desc.add("pdfSet", "CT14lo"s).setDescription("name of the PDF set to be used at the initial scale");
      desc.add("targetDIS", "isoscalar"s);
      return desc;
    }

  private:
    void eval() override {
      if (!q2_range_.contains(args_.q2) || args_.xbj < xbj_min_) {
        CG_WARNING("apfel:EvolutionStructureFunctions")
            << "(xBj=" << args_.xbj << ", Q^2=" << args_.q2 << ")"
            << " not in validity range (min.xBj = " << xbj_min_ << ", Q^2 = " << q2_range_ << ").";
        clear();
        return;
      }
      const auto q = std::sqrt(args_.q2);
      setF2(APFEL::StructureFunctionxQ(proc_, "F2", "total", args_.xbj, q));
      setFL(APFEL::StructureFunctionxQ(proc_, "FL", "total", args_.xbj, q));
    }

    const std::string proc_;
    const Limits q2_range_;
    const double xbj_min_;
  };
}  // namespace cepgen::apfel
using ApfelEvolutionStructureFunctions = cepgen::apfel::EvolutionStructureFunctions;
REGISTER_STRFUN("apfelEvol", 404, ApfelEvolutionStructureFunctions);
