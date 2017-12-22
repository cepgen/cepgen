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
        virtual double operator()( double q2, double xbj ) const = 0;

      protected:
        double xi( double q2, double xbj ) const;
    };

    class E143Ratio : public SigmaRatio
    {
      public:
        struct Parameterisation
        {
          double q2_b;
          std::array<double,6> a, b, c;
          static Parameterisation standard();
        };
        E143Ratio( const Parameterisation& param = Parameterisation::standard() ) : params_( param ) {}
        double operator()( double q2, double xbj ) const override;

      private:
        Parameterisation params_;
    };

    class CLASRatio : public SigmaRatio
    {
      public:
        CLASRatio() {}
        double operator()( double q2, double xbj ) const override;
    };
  }
}

#endif

