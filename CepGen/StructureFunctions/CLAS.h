#ifndef CepGen_StructureFunctions_CLAS_h
#define CepGen_StructureFunctions_CLAS_h

#include "StructureFunctions.h"
#include "CepGen/Physics/Constants.h"
#include <array>
#include <vector>

namespace CepGen
{
  namespace SF
  {
    class CLAS : public StructureFunctions
    {
      public:
        struct Parameterisation
        {
          static Parameterisation standard_neutron();
          static Parameterisation standard_proton();
          static Parameterisation standard_deuteron();

          enum { neutron = 0, proton = 1, deuteron = 2 } mode;
          double mp, mpi0;
          // SLAC fit parameters
          std::array<double,7> c_slac;
          // CLAS parameterisation
          double alpha, beta, dmu, dmup;
          std::array<double,3> x;
          std::array<double,4> b, ar, dmr, dgr;
          std::array<unsigned short,4> lr;
        };

        CLAS( const CLAS::Parameterisation& params = CLAS::Parameterisation::standard_proton() ) : params_( params ) {}

        CLAS operator()( double q2, double xbj ) const;

      private:
        void resbkg( double q2, double w, double& f2bkg, double& f2resn ) const;
        double f2slac( double q2, double xbj ) const;
        Parameterisation params_;
    };
  }
}

#endif
