/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#ifndef CepGen_FormFactors_Parameterisation_h
#define CepGen_FormFactors_Parameterisation_h

#include "CepGen/FormFactors/FormFactors.h"
#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Physics/HeavyIon.h"

namespace cepgen {
  class ParametersList;
  namespace strfun {
    class Parameterisation;
  }
  /// Form factors definition scope
  namespace formfac {
    /// Nucleon electromagnetic form factors parameterisation
    class Parameterisation : public NamedModule<std::string> {
    public:
      /// Steered parameterisation object constructor
      explicit Parameterisation(const ParametersList&);

      static ParametersDescription description();

      /// Dumping operator for standard output streams
      friend std::ostream& operator<<(std::ostream&, const Parameterisation*);
      /// Dumping operator for standard output streams
      friend std::ostream& operator<<(std::ostream&, const Parameterisation&);

      /// \f$\tau=Q^2/4m_p^2\f$ variable definition
      double tau(double q2) const;

      /// Compute all relevant form factors functions for a given \f$Q^2\f$ value
      virtual const FormFactors& operator()(double /*q2*/);

    protected:
      /// Proton magnetic moment
      static constexpr double MU = 2.79;

      /// Local form factors evaluation method
      virtual FormFactors compute(double) { return FormFactors{}; }

      const HeavyIon hi_;   ///< Incoming beam
      const double mass2_;  ///< Incoming beam squared mass
      const double mp_;     ///< Proton mass, in GeV/c\f$^2\f$
      const double mp2_;    ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$

    private:
      std::pair<double, FormFactors> last_value_{-1., FormFactors{}};
    };
    std::ostream& operator<<(std::ostream&, const FormFactors&);
  }  // namespace formfac
}  // namespace cepgen

#endif
