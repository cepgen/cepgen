/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2026  Laurent Forthomme
 *                2017-2019  Wolfgang Schaefer
 *                2019       Marta Luszczak
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

#ifndef CepGenProcesses_PPtoPair_h
#define CepGenProcesses_PPtoPair_h

#include "CepGen/Process/FactorisedProcess.h"

namespace cepgen::amplitudes {
  class HelicityAmplitudes;
}  // namespace cepgen::amplitudes

namespace cepgen::proc {
  /// Compute the 2-to-4 matrix element for a CE \f$\gamma\gamma\rightarrow Xf\bar X\f$ process
  class PPtoPair : public FactorisedProcess {
  public:
    explicit PPtoPair(const ParametersList& params, const std::vector<short>& hel_x1, const std::vector<short>& hel_x2);

    static ParametersDescription description();

  protected:
    void prepareFactorisedPhaseSpace() override;
    double offShellME() const;
    double helicityME() const;

    /// Parameters for the off-shell matrix element
    struct OffShellParameters : SteeredObject<OffShellParameters> {
      explicit OffShellParameters(const ParametersList&);

      static ParametersDescription description();

      int mat1{0}, mat2{0};
      int term_ll{0}, term_lt{0}, term_tt1{0}, term_tt2{0};
    } osp_;
    virtual double phiMatrixElement(
        double zp, double zm, double q2, const Momentum& vec_pho, const Momentum& vec_qt) const {
      (void)zp;
      (void)zm;
      (void)q2;
      (void)vec_pho;
      (void)vec_qt;
      return 0.;
    }

    static constexpr double kFourPi = 4. * M_PI;

    const std::vector<short> helicities_x1_;
    const std::vector<short> helicities_x2_;
    const enum class Mode { onShell = 0, offShell = 1, helicity = 2 } method_;
    std::unique_ptr<const amplitudes::HelicityAmplitudes> amplitudes_;  ///< helicity amplitudes (method 2 only)
    double mx2_{0.};
    double qx2_{0.};
    double colour_factor_{1.};
    std::function<double(double)> g_part1_{nullptr};
    std::function<double(double)> g_part2_{nullptr};

  private:
    static spdgids_t buildPdgIdsPair(pdgid_t pair) {
      return spdgids_t{static_cast<spdgid_t>(pair), -static_cast<spdgid_t>(pair)};
    }
    /// Coupling factor \f$N_c(4\pi\alpha_1Q^2)(4\pi\alpha_2Q^2)\f$ for given photon scales
    double couplingFactor(double q1, double q2) const;
  };
}  // namespace cepgen::proc
#endif
