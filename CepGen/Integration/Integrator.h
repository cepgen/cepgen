#ifndef CepGen_Integration_Integrator_h
#define CepGen_Integration_Integrator_h

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

#include <vector>
#include <memory>
#include <functional>

#include <string.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

namespace cepgen
{
  class Parameters;
  /// Monte-Carlo integrator instance
  class Integrator
  {
    public:
      /// Book the memory slots and structures for the integrator
      Integrator( const ParametersList& params );

      /**
       * Specify the function to be integrated
       * \param[in] ndim Number of dimensions on which the function will be integrated
       * \param[in] integrand Function to be integrated
       * \param[inout] params Run parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
       */
      void setFunction( unsigned int ndim, double integrand( double*, size_t, void* ), Parameters& params );
      /// Dimensional size of the phase space
      size_t size() const;
      /// Compute the function value at the given phase space point
      virtual double eval( const std::vector<double>& x );

      /**
       * Algorithm to perform the n-dimensional Monte Carlo integration of a given function.
       * \param[out] result_ The cross section as integrated for the given phase space restrictions
       * \param[out] abserr_ The uncertainty associated to the computed cross section
       */
      virtual void integrate( double& result_, double& abserr_ ) = 0;

      /// Algorithm name
      const std::string& name() const { return name_; }

      /// Random number generator instance
      const gsl_rng& rng() const { return *rng_; }
      /// Generate a uniformly distributed (between 0 and 1) random number
      double uniform() const;

    protected:
      const ParametersList params_;
      const std::string name_; ///< Integration algorithm name
      unsigned long seed_; ///< Random number generator seed
      struct gsl_rng_deleter
      {
        inline void operator()( gsl_rng* rng ) { gsl_rng_free( rng ); }
      };
      /// Instance of random number generator service
      std::unique_ptr<gsl_rng,gsl_rng_deleter> rng_;
      /// List of parameters to specify the integration range and the
      /// physics determining the phase space
      Parameters* input_params_;
      /// GSL structure storing the function to be integrated by this
      /// integrator instance (along with its parameters)
      std::unique_ptr<gsl_monte_function> function_;
      double result_, err_result_;
      bool initialised_;
  };
}

#endif
