/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

#ifndef CepGen_StructureFunctions_SuriYennie_h
#define CepGen_StructureFunctions_SuriYennie_h

#include "CepGen/StructureFunctions/Parameterisation.h"

namespace cepgen {
  namespace strfun {
    /// \f$F_{1,2,E,M}\f$ modelling by Suri and Yennie \cite Suri:1971yx
    class SuriYennie final : public Parameterisation {
    public:
      /// User-steered Suri-Yennie continuum structure functions calculator
      explicit SuriYennie(const ParametersList& params = ParametersList());
      static std::string description() { return "Suri-Yennie FE/FM"; }

      SuriYennie& eval(double xbj, double q2) override;

      double W1;  ///< Longitudinal form factor
      double W2;
      double FE;  ///< Electric proton form factor
      double FM;  ///< Magnetic proton form factor
    private:
      /// Collection of parameterisation-dependent couplings
      struct Parameters {
        double C1 = 0., C2 = 0.;
        double D1 = 0.;
        double rho2 = 0.;
        double Cp = 0., Bp = 0.;
      } sy_params_;
    };
  }  // namespace strfun
}  // namespace cepgen

#endif
