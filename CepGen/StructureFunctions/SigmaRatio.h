/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2017-2024  Laurent Forthomme
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

#ifndef CepGen_StructureFunctions_SigmaRatio_h
#define CepGen_StructureFunctions_SigmaRatio_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  /// A collector namespace for modellings of the \f$R=\sigma_L/\sigma_T\f$ ratio
  namespace sigrat {
    /// A generic modelling of the \f$R=\sigma_L/\sigma_T\f$ ratio
    class Parameterisation : public NamedModule<Parameterisation, std::string> {
    public:
      /// \f$R=\sigma_L/\sigma_T\f$ ratio computation algorithm constructor
      explicit Parameterisation(const ParametersList& params = ParametersList());
      /// Extract the longitudinal/transverse cross section ratio and associated error for a given \f$(x_{\rm Bj},Q^2)\f$ couple.
      virtual double operator()(double xbj, double q2, double& err) const = 0;

      static int index() { return 0; }  ///< Legacy integer index for sigma ratios modelling
      static ParametersDescription description();

    protected:
      /// \f$x_{\rm Bj}\f$ dependence for QCD-matching of R at high-\f$Q^2\f$
      static double theta(double xbj, double q2);
      const double mp_;   ///< Proton mass, in GeV/c\f$^2\f$
      const double mp2_;  ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$
    };
  }  // namespace sigrat
}  // namespace cepgen

#endif
