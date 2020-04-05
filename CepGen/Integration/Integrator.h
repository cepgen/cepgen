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
  class Integrand;
  /// Monte-Carlo integration algorithm
  class Integrator
  {
    public:
      /// Integrator algorithm constructor
      Integrator( const ParametersList& params );

      /// Specify the function to be integrated
      /// \param[in] integr Integrand object to be evaluated
      void setIntegrand( Integrand& integr );
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
      const int verbosity_; ///< Integrator verbosity
      /// A deleter object for GSL's random number generator
      struct gsl_rng_deleter
      {
        /// Destructor method for the random number generator service
        inline void operator()( gsl_rng* rng ) { gsl_rng_free( rng ); }
      };
      /// Instance of random number generator service
      std::unique_ptr<gsl_rng,gsl_rng_deleter> rng_;
      Integrand* integrand_; ///< Integrand to be evaluated
      /// GSL structure storing the function to be integrated by this
      /// integrator instance (along with its parameters)
      std::unique_ptr<gsl_monte_function> function_;
      double result_; ///< Result of the last integration
      double err_result_; ///< Standard deviation for the last integration
      bool initialised_; ///< Has the algorithm alreay been initialised?

    private:
      /**
       * \brief A GSL wrapper to define a functor as an integrable functional
       * \tparam F Class member signature
       */
      template<typename F>
      class gsl_monte_function_wrapper : public gsl_monte_function
      {
        public:
          gsl_monte_function_wrapper( const F& func, size_t ndim ) :
            func_( func ) {
            f = &gsl_monte_function_wrapper::eval;
            dim = ndim;
            params = this;
          }

        private:
          /// Static integrable functional
          static double eval( double* x, size_t ndim, void* params ) {
            return static_cast<gsl_monte_function_wrapper*>( params )->func_( x, ndim, params );
          }
          /// Reference to the functor
          const F& func_;
      };
  };
}

#endif
