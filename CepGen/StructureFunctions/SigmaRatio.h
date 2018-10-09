#ifndef CepGen_StructureFunctions_SigmaRatio_h
#define CepGen_StructureFunctions_SigmaRatio_h

#include <array>

namespace CepGen
{
  namespace SF
  {
    /// A generic modelling of the \f$R=\sigma_L/\sigma_T\f$ ratio
    class SigmaRatio
    {
      public:
        SigmaRatio() {}
        /// Extract the longitudinal/transverse cross section ratio and associated error for a given \f$(x_{\rm Bj},Q^2)\f$ couple.
        virtual double operator()( double xbj, double q2, double& err ) const = 0;

      protected:
        /// \f$x_{\rm Bj}\f$ dependence for QCD-matching of R at high-\f$Q^2\f$
        double theta( double xbj, double q2 ) const;
        static const double mp_; ///< Proton mass, in GeV/c\f${}^2\f$
        static const double mp2_; ///< Squared proton mass, in GeV\f${}^2\f$/c\f${}^4\f$
    };

    /// E143 experimental R measurement \cite Abe:1998ym
    class E143Ratio : public SigmaRatio
    {
      public:
        struct Parameterisation
        {
          double q2_b, lambda2;
          std::array<double,6> a, b, c;
          static Parameterisation standard();
        };
        explicit E143Ratio( const Parameterisation& param = Parameterisation::standard() );
        double operator()( double xbj, double q2, double& err ) const override;

      private:
        Parameterisation params_;
    };

    /** \brief SLAC experimental R measurement \cite Whitlow:1990gk
     * \warning valid for \f$Q^2\f$ > 0.3 GeV\f${}^2\f$
     */
    class R1990Ratio: public SigmaRatio
    {
      public:
        struct Parameterisation
        {
          double lambda2;
          std::array<double,3> b;
          static Parameterisation standard();
        };
        explicit R1990Ratio( const Parameterisation& param = Parameterisation::standard() );
        double operator()( double xbj, double q2, double& err ) const override;

      private:
        Parameterisation params_;
    };

    /// CLAS experimental R measurement
    class CLASRatio : public SigmaRatio
    {
      public:
        CLASRatio() {}
        double operator()( double xbj, double q2, double& err ) const override;
    };

    /// Sibirtsev & Blunden parameterisation of the R ratio \cite Sibirtsev:2013cga
    class SBRatio : public SigmaRatio
    {
      public:
        SBRatio() {}
        double operator()( double xbj, double q2, double& err ) const override;
    };
  }
}

#endif
