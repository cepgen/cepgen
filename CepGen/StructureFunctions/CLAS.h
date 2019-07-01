#ifndef CepGen_StructureFunctions_CLAS_h
#define CepGen_StructureFunctions_CLAS_h

#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <array>
#include <vector>

namespace cepgen
{
  namespace strfun
  {
    /// \brief CLAS parameterisation for nucleon data at \f$Q^2\f$ > 0.5 GeV\f$^2\f$ and \f$x_{\rm Bj}\f$ > 0.15
    /// \note This code was provided on 2016-04-13 by Silvano Simula and reflects the parameterisation used in \cite Osipenko:2003bu (CLAS) and described in \cite Ricco:1998yr.
    class CLAS : public Parameterisation
    {
      public:
        /// List of steering parameters for a physics case
        struct Parameters
        {
          /// Standard parameterisation of a parton-from-neutron emission
          static Parameters standard_neutron();
          /// Standard parameterisation of a parton-from-proton emission
          static Parameters standard_proton();
          /// Standard parameterisation of a parton-from-deuteron emission
          static Parameters standard_deuteron();

          /// Physical properties associated to a resonance
          struct Resonance
          {
            double amplitude, mass, width;
            short angular_momentum;
          };

          enum { neutron = 0, proton = 1, deuteron = 2 } mode; ///< Nucleon type
          double mp; ///< Proton mass
          double mpi0; ///< Neutral pion mass
          // SLAC fit parameters
          std::array<double,7> c_slac;
          // CLAS parameterisation
          double alpha, beta, mu, mup;
          std::array<double,3> x;
          std::array<double,4> b;
          std::vector<Resonance> resonances;
          std::array<unsigned short,4> lr;
        };

        /// Standard parameterisation interpolator constructor (photon from proton)
        explicit CLAS( const ParametersList& params = ParametersList() );

        CLAS& operator()( double xbj, double q2 ) override;

      private:
        /// \brief Method to evaluate the background/resonance terms of
        ///  the modulating function for the nucleon
        /// \note SLAC parameterisation
        std::pair<double,double> resbkg( double q2, double w ) const;
        /// \brief Method to evaluate the deep inelastic structure function
        /// \f$F_{2}^{N}\f$ using the SLAC parameterisation
        /// \param[in] q2 squared four-momentum transfer in GeV\f$^2\f$
        /// \param[in] xbj Bjorken scaling variable
        /// \return \f$F_{2}^{N}\f$
        double f2slac( double xbj, double q2 ) const;
        Parameters params_;
        static constexpr double COEFF = 6.08974;
    };
  }
}

#endif
