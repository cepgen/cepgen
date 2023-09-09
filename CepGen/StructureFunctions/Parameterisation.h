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

#ifndef CepGen_StructureFunctions_Parameterisation_h
#define CepGen_StructureFunctions_Parameterisation_h

#include <iosfwd>
#include <memory>

#include "CepGen/Modules/NamedModule.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

namespace cepgen {
  /// Structure functions modelling scope
  namespace strfun {
    /// Base object for the parameterisation of nucleon structure functions
    class Parameterisation : public NamedModule<int> {
    public:
      /// User-steered parameterisation object constructor
      explicit Parameterisation(const ParametersList&);
      virtual ~Parameterisation() = default;

      /// Generic description for the structure functions
      static ParametersDescription description();

      /// Human-readable dump of the SF parameterisation at this (xBj,Q^2) value
      friend std::ostream& operator<<(std::ostream&, const Parameterisation&);

      /// Longitudinal/transverse cross section ratio parameterisation used to compute \f$F_{1/L}\f$
      const sigrat::Parameterisation* sigmaRatio() const { return r_ratio_.get(); }

      /// Compute all relevant structure functions for a given \f$(x_{\rm Bj},Q^2)\f$ couple
      /// \param[in] xbj Bjorken's x variable
      /// \param[in] q2 Squared 4-momentum transfer (in GeV^2)
      Parameterisation& operator()(double /*xbj*/, double /*q2*/);
      /// Transverse structure function
      double F2(double xbj, double q2);
      /// Longitudinal structure function
      double FL(double xbj, double q2);
      /// Longitudinal form factor
      double W1(double xbj, double q2);
      double W2(double xbj, double q2);
      /// Electric proton form factor
      double FE(double xbj, double q2);
      /// Magnetic proton form factor
      double FM(double xbj, double q2);
      /// \f$F_1\f$ structure function
      double F1(double xbj, double q2);

      struct Arguments {
        bool operator==(const Arguments& oth) const { return xbj == oth.xbj && q2 == oth.q2; }
        bool valid() const { return q2 >= 0. && xbj >= 0. && xbj < 1.; }
        friend std::ostream& operator<<(std::ostream&, const Arguments&);
        double xbj{-1.}, q2{-1.};
      };
      struct Values {
        void clear() { f2 = 0., fl = 0., w1 = 0., w2 = 0., f2 = 0., fm = 0.; }
        friend std::ostream& operator<<(std::ostream&, const Values&);
        double f2{0.};  ///< Last computed transverse structure function value
        double fl{0.};  ///< Last computed longitudinal structure function value
        // alternative quantities
        double w1{0.};  ///< Longitudinal form factor
        double w2{0.};
        double fe{0.};  ///< Electric proton form factor
        double fm{0.};  ///< Magnetic proton form factor
      };

    protected:
      /// Local structure functions evaluation method
      virtual void eval() = 0;
      /// Compute the longitudinal structure function for a given point
      virtual Parameterisation& computeFL(double xbj, double q2);
      /// Compute the longitudinal structure function for a given point
      virtual Parameterisation& computeFL(double xbj, double q2, double r);

      /// Reset the structure functions values
      Parameterisation& clear();
      //-- fill in the structure functions values
      Parameterisation& setF1F2(double f1, double f2);
      Parameterisation& setF2(double f2);
      Parameterisation& setFL(double fl);
      Parameterisation& setW1(double w1);
      Parameterisation& setW2(double w2);
      Parameterisation& setFE(double fe);
      Parameterisation& setFM(double fm);

      /// Compute the dimensionless variable \f$\tau=\frac{4x_{\rm Bj}^2m_p^2}{Q^2}\f$
      double tau(double xbj, double q2) const;
      /// Dimensionless variable \f$\gamma^2=1+\frac{4x_{\rm Bj}^m_p^2}{Q^2}=1+\tau\f$
      double gamma2(double xbj, double q2) const;

    private:
      /// Longitudinal/transverse cross section ratio parameterisation used to compute \f$F_{1/L}\f$
      const std::unique_ptr<sigrat::Parameterisation> r_ratio_;

    protected:
      const double mp_;      ///< Proton mass, in GeV/c^2
      const double mp2_;     ///< Squared proton mass, in GeV^2/c^4
      const double mx_min_;  ///< Minimum diffractive mass, in GeV/c^2

      Arguments args_;  ///< Last \f$(x_{\rm Bj},Q^2)\f$ couple computed

    private:
      Values vals_;
      bool fl_computed_{false};
    };
  }  // namespace strfun
}  // namespace cepgen

#endif
