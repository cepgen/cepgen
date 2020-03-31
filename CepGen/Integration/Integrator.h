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
  /// Monte-Carlo integration algorithm
  class Integrator
  {
    public:
      /// Integrator algorithm constructor
      Integrator( const ParametersList& params );

      /// Specify the function to be integrated
      /// \param[in] ndim Number of dimensions on which the function will be integrated
      /// \param[in] integrand Function to be integrated
      /// \param[inout] params Run parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
      void setFunction( unsigned int ndim, double integrand( double*, size_t, void* ), Parameters& params );
      /// Dimensional size of the phase space
      size_t size() const;
      /// Compute the function value at the given phase space point
      virtual double eval( const std::vector<double>& x );

      /// Perform the multidimensional Monte Carlo integration
      /// \param[out] result_ The cross section as integrated for the given phase space restrictions
      /// \param[out] abserr_ The uncertainty associated to the computed cross section
      virtual void integrate( double& result_, double& abserr_ ) = 0;

      /// Integration algorithm name
      const std::string& name() const { return name_; }

      /// Random number generator instance
      const gsl_rng& rng() const { return *rng_; }
      /// Generate a uniformly distributed (between 0 and 1) random number
      double uniform() const;

    protected:
      const ParametersList params_; ///< Steering parameters for this algorithm
      const std::string name_; ///< Integration algorithm name
      unsigned long seed_; ///< Random number generator seed
      /// A deleter object for GSL's random number generator
      struct gsl_rng_deleter
      {
        /// Destructor method for the random number generator service
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
      double result_; ///< Result of the last integration
      double err_result_; ///< Standard deviation for the last integration
      bool initialised_; ///< Has the algorithm alreay been initialised?
  };
}

#endif
