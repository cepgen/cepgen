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

#include <APFEL/APFEL.h>

#include <cmath>

#include "CepGen/Core/Exception.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

using namespace std::string_literals;

namespace cepgen {
  namespace apfel {
    class EvolutionStructureFunctions final : public strfun::Parameterisation {
    public:
      explicit EvolutionStructureFunctions(const ParametersList& params)
          : strfun::Parameterisation(params),
            proc_(steer<std::string>("processDIS")),
            xbj_min_(steer<double>("xBjmin")) {
        const auto q2range = steer<Limits>("q2range");
        if (!q2range.valid())
          throw CG_FATAL("apfel:EvolutionStructureFunctions") << "Invalid Q^2 range: " << q2range << ".";
        const auto qrange = q2range.compute([](double lim) { return std::sqrt(lim); });

        APFEL::SetMassScheme(steer<std::string>("massScheme"));
        APFEL::SetProcessDIS(proc_);
        APFEL::SetQLimits(qrange.min(), qrange.max());
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
        desc.add("q2range", Limits{1., 1.e6});
        desc.add<double>("xBjmin", 1.e-6);
        desc.add("massScheme", "ZM-VFNS"s);
        desc.add("processDIS", "NC"s);
        desc.add<int>("maxFlavourAlpha", 5);
        desc.add<int>("maxFlavourPDFs", 5);
        desc.add("pdfSet", "CT14nnlo"s);
        desc.add("targetDIS", "isoscalar"s);
        return desc;
      }

    private:
      void eval() override {
        if (args_.xbj < xbj_min_) {
          clear();
          return;
        }
        const auto q = std::sqrt(args_.q2);
        setF2(APFEL::StructureFunctionxQ(proc_, "F2", "total", args_.xbj, q));
        setFL(APFEL::StructureFunctionxQ(proc_, "FL", "total", args_.xbj, q));
      }

      const std::string proc_;
      const double xbj_min_;
    };
  }  // namespace apfel
}  // namespace cepgen
using ApfelEvolutionStructureFunctions = cepgen::apfel::EvolutionStructureFunctions;
REGISTER_STRFUN("apfelEvol", 404, ApfelEvolutionStructureFunctions);
