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

#ifndef CepGen_StructureFunctions_Parameterisation_h
#define CepGen_StructureFunctions_Parameterisation_h

#include <iosfwd>
#include <memory>

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  namespace sigrat {
    class Parameterisation;
  }
  /// Structure functions modelling scope
  namespace strfun {
    /// Proton structure function to be used in the outgoing state description
    /// \note Values correspond to the LPAIR legacy steering card values
    enum struct Type {
      Invalid = 0,
      Electron = 1,
      ElasticProton = 2,
      SuriYennie = 11,
      SzczurekUleshchenko = 12,
      BlockDurandHa = 13,
      SuriYennieAlt = 14,
      FioreBrasse = 101,
      ChristyBosted = 102,
      CLAS = 103,
      FioreBrasseAlt = 104,
      ALLM91 = 201,
      ALLM97 = 202,
      GD07p = 203,
      GD11p = 204,
      MSTWgrid = 205,
      HHT_ALLM = 206,
      HHT_ALLM_FT = 207,
      Schaefer = 301,
      Shamov = 302,
      KulaginBarinov = 303,
      Partonic = 401,
    };
    /// Human-readable description of this SF parameterisation type
    std::ostream& operator<<(std::ostream&, const strfun::Type&);
    /// Generic placeholder for the parameterisation of nucleon structure functions
    class Parameterisation : public NamedModule<int> {
    public:
      /// Standard SF parameterisation constructor
      explicit Parameterisation(double f2 = 0., double fl = 0.);
      /// Copy constructor
      Parameterisation(const Parameterisation&);
      /// User-steered parameterisation object constructor
      explicit Parameterisation(const ParametersList&);
      virtual ~Parameterisation() = default;

      /// Assign from another SF parameterisation object
      Parameterisation& operator=(const Parameterisation& sf);

      /// Generic description for the structure functions
      static ParametersDescription description();

      /// Human-readable dump of the SF parameterisation at this (xBj,Q^2) value
      friend std::ostream& operator<<(std::ostream&, const Parameterisation*);
      /// Human-readable dump of the SF parameterisation at this (xBj,Q^2) value
      friend std::ostream& operator<<(std::ostream&, const Parameterisation&);
      /// Human-readable description of this SF parameterisation
      virtual std::string describe() const;  ///< Human-readable description of this SF set

      /// Longitudinal/transverse cross section ratio parameterisation used to compute \f$F_{1/L}\f$
      const sigrat::Parameterisation* sigmaRatio() const { return r_ratio_.get(); }

      /// Compute all relevant structure functions for a given \f$(x_{\rm Bj},Q^2)\f$ couple
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
        bool valid() const { return q2 >= 0. && xbj >= 0. && xbj <= 1.; }
        friend std::ostream& operator<<(std::ostream&, const Arguments&);
        double xbj{-1.}, q2{-1.};
      };

    protected:
      /// Local structure functions evaluation method
      /// \param[in] xbj Bjorken's x variable
      /// \param[in] q2 Squared 4-momentum transfer (in GeV^2)
      virtual Parameterisation& eval(double xbj, double q2);
      /// Compute the longitudinal structure function for a given point
      virtual Parameterisation& computeFL(double xbj, double q2);
      /// Compute the longitudinal structure function for a given point
      virtual Parameterisation& computeFL(double xbj, double q2, double r);

      //-- fill in the structure functions values
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

      const double mp_;      ///< Proton mass, in GeV/c^2
      const double mp2_;     ///< Squared proton mass, in GeV^2/c^4
      const double mx_min_;  ///< Minimum diffractive mass, in GeV/c^2

    private:
      Arguments old_vals_;  ///< Last \f$(x_{\rm Bj},Q^2)\f$ couple computed
      /// Longitudinal/transverse cross section ratio parameterisation used to compute \f$F_{1/L}\f$
      std::shared_ptr<sigrat::Parameterisation> r_ratio_;
      double f2_{0.};  ///< Last computed transverse structure function value
      double fl_{0.};  ///< Last computed longitudinal structure function value
      bool fl_computed_{false};
      // alternative quantities
      double w1_{0.};  ///< Longitudinal form factor
      double w2_{0.};
      double fe_{0.};  ///< Electric proton form factor
      double fm_{0.};  ///< Magnetic proton form factor
    };
  }  // namespace strfun
}  // namespace cepgen

#endif
