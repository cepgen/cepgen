/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  namespace formfac {
    /// \cite A1:2013fsc
    class A1Elastic final : public Parameterisation {
    public:
      explicit A1Elastic(const ParametersList& params)
          : Parameterisation(params),
            coeff_e_(steer<std::vector<double> >("coeffE")),
            coeff_m_(steer<std::vector<double> >("coeffM")),
            max_interp_q2_(steer<double>("q2interp")) {
        if (coeff_e_.size() < 3)
          throw CG_FATAL("A1Elastic") << "Invalid coefficients multiplicity for the G_E functional form!";
        if (coeff_m_.size() < 3)
          throw CG_FATAL("A1Elastic") << "Invalid coefficients multiplicity for the G_M functional form!";
        const auto grid_filename = steerPath("A1SplinesGrid");
        std::ifstream grid_file(grid_filename);
        for (std::string line; std::getline(grid_file, line);) {
          const auto vals = utils::split(line, ' ');
          if (vals.size() < 6)  // should be 13
            continue;
          const auto q2 = std::stod(vals.at(0)), ge = std::stod(vals.at(1)), gm = std::stod(vals.at(5));
          coh_grid_.insert({q2}, {ge, gm});
        }
        coh_grid_.init();
        min_interp_q2_ = coh_grid_.min().at(0);
        CG_DEBUG("A1Elastic") << "Splines interpolation grid file loaded from '" << grid_filename << ". "
                              << "Q^2 range: " << coh_grid_.boundaries().at(0) << " GeV^2.";
      }

      static ParametersDescription description() {
        auto desc = Parameterisation::description();
        desc.setDescription("A1 elastic");
        desc.add<std::vector<double> >("coeffE", {0.98462, 0.68414, 0.01933})
            .setDescription("coefficients for the G_E functional form");
        desc.add<std::vector<double> >("coeffM", {0.28231, 1.34919, 0.55473})
            .setDescription("coefficients for the G_M functional form");
        desc.add<std::string>("A1SplinesGrid", "External/PhysRevC.90.015206.SplinesWithVariableKnots.dat");
        desc.add<double>("q2interp", 10.).setDescription("maximal Q^2 at which interpolation is performed (in GeV^2)");
        return desc;
      }

    private:
      /// Friedrich-Walcher double dipole parameterisation
      /// \cite Friedrich:2003iz
      double doubleDipoleGEM(double q2, const std::vector<double>& coeffs) const {
        return coeffs.at(0) * std::pow(1. + q2 / coeffs.at(1), -2) +
               (1. - coeffs.at(0)) * std::pow(1. + q2 / coeffs.at(2), -2);
      }
      FormFactors compute(double q2) override {
        FormFactors out;
        if (q2 < min_interp_q2_) {
          const auto& min_vals = coh_grid_.values().begin()->second;
          out.GE = (1. + q2 * (min_vals.at(0) - 1.) / 0.005);
          out.GM = (1. + q2 * (min_vals.at(1) - 1.) / 0.005) * MU;
          return out;
        }
        if (q2 < max_interp_q2_) {
          const auto& grid_vals = coh_grid_.eval({q2});
          out.GE = grid_vals.at(0);
          out.GM = grid_vals.at(1) * MU;
          return out;
        }
        out.GE = doubleDipoleGEM(q2, coeff_e_);
        out.GM = doubleDipoleGEM(q2, coeff_m_) * MU;
        return out;
      }
      const std::vector<double> coeff_e_, coeff_m_;
      GridHandler<1, 2> coh_grid_{GridType::linear};
      double min_interp_q2_{0.};
      const double max_interp_q2_;
    };
  }  // namespace formfac
}  // namespace cepgen
using cepgen::formfac::A1Elastic;
REGISTER_FORMFACTORS("A1Elastic", A1Elastic);
