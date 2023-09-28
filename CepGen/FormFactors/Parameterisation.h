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
      friend std::ostream& operator<<(std::ostream&, const Parameterisation&);

      /// \f$\tau=Q^2/4m_p^2\f$ variable definition
      double tau(double q2) const;

      /// Compute all relevant form factors functions for a given \f$Q^2\f$ value
      virtual const FormFactors& operator()(double /*q2*/);

    protected:
      /// Proton magnetic moment
      static constexpr double MU = 2.792847337;

      /// Local form factors evaluation method
      virtual void eval() = 0;
      /// Set the form factors directly
      void setFEFM(double fe, double fm);
      /// Set the Sachs form factors
      void setGEGM(double ge, double gm);

      const pdgid_t pdg_id_;  ///< Incoming beam
      const double mass2_;    ///< Incoming beam squared mass
      const double mp_;       ///< Proton mass, in GeV/c\f$^2\f$
      const double mp2_;      ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$

      /// Virtuality at which the form factors are evaluated
      double q2_{-1.};
      /// Last form factors computed
      FormFactors ff_{};
    };
    std::ostream& operator<<(std::ostream&, const FormFactors&);
  }  // namespace formfac
}  // namespace cepgen

#endif
