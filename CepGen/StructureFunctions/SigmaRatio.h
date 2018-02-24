#ifndef CepGen_StructureFunctions_SigmaRatio_h
#define CepGen_StructureFunctions_SigmaRatio_h

#include <array>

namespace CepGen
{
  namespace SF
  {
    class SigmaRatio
    {
      public:
        SigmaRatio() {}
        /// Extract the longitudinal/transverse cross section ratio and associated error for a given Q²/\f$x_{\textrm{Bj}}\f$ couple.
        virtual double operator()( double q2, double xbj, double& err ) const = 0;

      protected:
        /// \f$x_{\textrm{Bj}}\f$ dependence for QCD-matching of R at high-Q²
        double theta( double q2, double xbj ) const;
    };

    // Reference: arXiv:hep-ex/9808028
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
        double operator()( double q2, double xbj, double& err ) const override;

      private:
        Parameterisation params_;
    };

    /// \warning valid for Q² > 0.3 GeV²
    // Reference: Phys.Lett. B 250 (1990) 193-198 (https://inspirehep.net/record/296980)
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
        double operator()( double q2, double xbj, double& err ) const override;

      private:
        Parameterisation params_;
    };

    class CLASRatio : public SigmaRatio
    {
      public:
        CLASRatio() {}
        double operator()( double q2, double xbj, double& err ) const override;
    };

    /// Sibirtsev & Blunden parameterisation of the R ratio
    // Reference: Phys.Rev. C 88,065202 (2013)
    class SBRatio : public SigmaRatio
    {
      public:
        SBRatio() {}
        double operator()( double q2, double xbj, double& err ) const override;
    };
  }
}

#endif

