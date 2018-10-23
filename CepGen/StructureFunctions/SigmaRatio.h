#ifndef CepGen_StructureFunctions_SigmaRatio_h
#define CepGen_StructureFunctions_SigmaRatio_h

#include "CepGen/Core/ParametersList.h"

#include <memory>

namespace cepgen
{
  /// A collector namespace for modellings of the \f$R=\sigma_L/\sigma_T\f$ ratio
  namespace sigrat
  {
    enum struct Type { E143 = 1, R1990 = 2, CLAS = 3, SibirtsevBlunden = 4 };
    /// A generic modelling of the \f$R=\sigma_L/\sigma_T\f$ ratio
    class Parameterisation
    {
      public:
        Parameterisation( const ParametersList& params = ParametersList() );
        /// Extract the longitudinal/transverse cross section ratio and associated error for a given \f$(x_{\rm Bj},Q^2)\f$ couple.
        virtual double operator()( double xbj, double q2, double& err ) const = 0;

        /// Build a R-ratio parameterisation for a given type
        static std::shared_ptr<Parameterisation> build( const ParametersList& );

      protected:
        /// \f$x_{\rm Bj}\f$ dependence for QCD-matching of R at high-\f$Q^2\f$
        double theta( double xbj, double q2 ) const;
        static const double mp_; ///< Proton mass, in GeV/c\f$^2\f$
        static const double mp2_; ///< Squared proton mass, in GeV\f$^2\f$/c\f$^4\f$
    };

    /// E143 experimental R measurement \cite Abe:1998ym
    class E143 : public Parameterisation
    {
      public:
        explicit E143( const ParametersList& params = ParametersList() );
        double operator()( double xbj, double q2, double& err ) const override;

      private:
        double q2_b_, lambda2_;
        std::vector<double> a_, b_, c_;
    };

    /** \brief SLAC experimental R measurement \cite Whitlow:1990gk
     * \warning valid for \f$Q^2\f$ > 0.3 GeV\f$^2\f$
     */
    class R1990: public Parameterisation
    {
      public:
        explicit R1990( const ParametersList& params = ParametersList() );
        double operator()( double xbj, double q2, double& err ) const override;

      private:
        double lambda2_;
        std::vector<double> b_;
    };

    /// CLAS experimental R measurement
    class CLAS : public Parameterisation
    {
      public:
        explicit CLAS( const ParametersList& params = ParametersList() );
        double operator()( double xbj, double q2, double& err ) const override;

      private:
        std::vector<double> p_;
        double wth_, q20_;
    };

    /// Sibirtsev & Blunden parameterisation of the R ratio \cite Sibirtsev:2013cga
    class SibirtsevBlunden : public Parameterisation
    {
      public:
        explicit SibirtsevBlunden( const ParametersList& params = ParametersList() );
        double operator()( double xbj, double q2, double& err ) const override;

      private:
        double a_, b1_, b2_, c_;
    };
  }
}

#endif
