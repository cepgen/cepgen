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

#ifndef CepGen_FormFactors_Parameterisation_h
#define CepGen_FormFactors_Parameterisation_h

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/Physics/Modes.h"

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
      /// Empty parameterisation object constructor
      explicit Parameterisation();
      /// Steered parameterisation object constructor
      explicit Parameterisation(const ParametersList&);
      /// Copy constructor
      Parameterisation(const Parameterisation&);

      static ParametersDescription description();

      /// Dumping operator for standard output streams
      friend std::ostream& operator<<(std::ostream&, const Parameterisation*);
      /// Dumping operator for standard output streams
      friend std::ostream& operator<<(std::ostream&, const Parameterisation&);

      /// Specify the structure functions modelling where applicable
      void setStructureFunctions(strfun::Parameterisation* sfmod) { str_fun_ = sfmod; }
      /// Retrieve the structure function interpolation object used here
      strfun::Parameterisation* structureFunctions() const { return str_fun_; }

      /// \f$\tau=Q^2/4m_p^2\f$ variable definition
      double tau(double q2) const;

      /// Compute all relevant form factors functions for a given \f$Q^2\f$ value
      Parameterisation& operator()(const mode::Beam& /*type*/, double /*q2*/, double mf2 = 0.);

    protected:
      /// Proton magnetic moment
      static constexpr double MU = 2.79;

      /// Local form factors evaluation method
      virtual void compute(double) {}

      const double mp_;   ///< Proton mass, in GeV/c\f$^2\f$
      const double mp2_;  ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$

    private:
      strfun::Parameterisation* str_fun_{nullptr};
      double last_q2_{-1.};

    public:
      double FE{0.};  ///< Electric form factor
      double FM{0.};  ///< Magnetic form factor

      double GE{0.};  ///< Sachs electric form factor
      double GM{0.};  ///< Sachs magnetic form factor
    };
  }  // namespace formfac
}  // namespace cepgen

#endif
