/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

namespace cepgen {
  class ParametersList;
  class NachtmannAmplitudes {
  public:
    enum class Mode { SM, W, Wbar, phiW, WB };
    NachtmannAmplitudes(const ParametersList&);
    double operator()(double shat, double that, double uhat, short lam1, short lam2, short lam3, short lam4) const;
    const Mode& mode() const { return mode_; }

  private:
    Mode mode_;
    struct EFTParameters {
      explicit EFTParameters(const ParametersList& params);
      const double s1, mH;
      double c1() const { return sqrt(1. - s1 * s1); }
    } eft_ext_;
    const double mw2_;
  };
}  // namespace cepgen
