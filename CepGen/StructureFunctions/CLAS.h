#ifndef CepGen_StructureFunctions_CLAS_h
#define CepGen_StructureFunctions_CLAS_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <array>
#include <vector>

namespace CepGen
{
  namespace SF
  {
    /// CLAS parameterisation developed to describe nucleon data at \f$Q^2 > 0.5\f$ GeV² and \f$x_{Bj} > 0.15\f$.
    /// \note This code was provided on 2016-04-13 by Silvano Simula and
    ///  reflects the parametrisation used in hep-ph/0301204 (CLAS) and
    ///  described in hep-ph/9901360.
    class CLAS : public StructureFunctions
    {
      public:
        struct Parameterisation
        {
          static Parameterisation standard_neutron();
          static Parameterisation standard_proton();
          static Parameterisation standard_deuteron();

          /// Physical properties associated to a resonance
          struct Resonance
          {
            double amplitude, mass, width;
            short angular_momentum;
          };

          enum { neutron = 0, proton = 1, deuteron = 2 } mode;
          double mp, mpi0;
          // SLAC fit parameters
          std::array<double,7> c_slac;
          // CLAS parameterisation
          double alpha, beta, mu, mup;
          std::array<double,3> x;
          std::array<double,4> b;
          std::vector<Resonance> resonances;
          std::array<unsigned short,4> lr;
        };

        explicit CLAS( const CLAS::Parameterisation& params = CLAS::Parameterisation::standard_proton() );

        CLAS& operator()( double q2, double xbj ) override;

      private:
        /// \brief Method to evaluate the background/resonance terms of
        ///  the modulating function for the nucleon
        /// \note SLAC parameterisation
        std::pair<double,double> resbkg( double q2, double w ) const;
        /// \brief Method to evaluate the deep inelastic structure function
        /// \f$F_{2}^{N}\f$ using the SLAC parameterisation
        /// \param[in] q2 squared four-momentum transfer in GeV²
        /// \param[in] xbj Bjorken scaling variable
        /// \return \f$F_{2}^{N}\f$
        double f2slac( double q2, double xbj ) const;
        Parameterisation params_;
    };
  }
}

#endif
