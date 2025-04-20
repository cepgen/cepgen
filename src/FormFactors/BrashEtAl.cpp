/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Modules/FormFactorsFactory.h"

using namespace std::string_literals;

namespace cepgen::formfac {
  /// \cite Brash:2001qq
  class BrashEtAl final : public Parameterisation {
  public:
    explicit BrashEtAl(const ParametersList& params)
        : Parameterisation(params),
          gm_coefficient_(steer<std::vector<double> >("coeffGM"s)),
          r_coefficient_(steer<std::vector<double> >("coeffR"s)),
          max_q2_(steer<double>("q2max")) {
      if (gm_coefficient_.size() != 5)
        throw CG_FATAL("BrashEtAl") << "Invalid coefficients multiplicity for the G_M functional form!";
      if (r_coefficient_.size() != 2)
        throw CG_FATAL("BrashEtAl") << "Invalid coefficients multiplicity for the G_E/G_M ratio functional form!";
    }

    static ParametersDescription description();

  private:
    void eval() override {
      if (q2_ > max_q2_)
        CG_WARNING("BrashEtAl") << "Q² = " << q2_ << " > " << max_q2_ << " GeV² = max(Q²).\n\t"
                                << "Brash et al. FF parameterisation not designed for high-Q² values.";
      const auto r = std::min(1., 1. - r_coefficient_.at(0) * (q2_ - r_coefficient_.at(1)));
      if (r < 0.)
        return;
      const auto q = std::sqrt(q2_),
                 gm = 1. / (1. + q * (gm_coefficient_.at(0) +
                                      q * (gm_coefficient_.at(1) +
                                           q * (gm_coefficient_.at(2) +
                                                q * (gm_coefficient_.at(3) + q * gm_coefficient_.at(4))))));

      setGEGM(r * gm, MU * gm);
    }
    const std::vector<double> gm_coefficient_, r_coefficient_;
    const double max_q2_;
  };

  ParametersDescription BrashEtAl::description() {
    auto desc = Parameterisation::description();
    desc.setDescription("Brash et al.");
    desc.add("coeffGM"s, std::vector{0.116, 2.874, 0.241, 1.006, 0.345})
        .setDescription("coefficients for the G_M functional form");
    desc.add("coeffR"s, std::vector{0.13, 0.04}).setDescription("coefficients for the G_E/G_M ratio functional form");
    desc.add("q2max"s, 7.7).setDescription("maximal Q^2 supported (in GeV^2)");
    return desc;
  }
}  // namespace cepgen::formfac
using cepgen::formfac::BrashEtAl;
REGISTER_FORMFACTORS("Brash", BrashEtAl);
