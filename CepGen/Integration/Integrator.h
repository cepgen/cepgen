#ifndef CepGen_Integration_Integrator_h
#define CepGen_Integration_Integrator_h

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Event/Event.h"

#include <vector>
#include <random>

#include <string.h>

namespace cepgen
{
  class Integrand;
  /// Monte-Carlo integration algorithm
  class Integrator
  {
    public:
      /// Integrator algorithm constructor
      Integrator( const ParametersList& params );

      /// Specify the function to be integrated
      /// \param[in] integr Integrand object to be evaluated
      virtual void setIntegrand( Integrand& integr );
      /// Dimensional size of the phase space
      size_t size() const;
      /// Integration algorithm name
      const std::string& name() const { return name_; }

      /// Compute the function value at the given phase space point
      virtual double eval( const std::vector<double>& x ) const;
      /// Generate a uniformly distributed (between 0 and 1) random number
      virtual double uniform() const;

      /// Perform the multidimensional Monte Carlo integration
      /// \param[out] result_ The cross section as integrated for the given phase space restrictions
      /// \param[out] abserr_ The uncertainty associated to the computed cross section
      virtual void integrate( double& result_, double& abserr_ ) = 0;

    protected:
      const ParametersList params_; ///< Steering parameters for this algorithm
      const std::string name_; ///< Integration algorithm name
      const unsigned long seed_; ///< Random number generator seed
      int verbosity_; ///< Integrator verbosity
      Integrand* integrand_; ///< Integrand to be evaluated
      double result_; ///< Result of the last integration
      double err_result_; ///< Standard deviation for the last integration
      bool initialised_; ///< Has the algorithm alreay been initialised?
  };
}

#endif
