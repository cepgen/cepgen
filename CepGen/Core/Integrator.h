#ifndef CepGen_Core_Integrator_h
#define CepGen_Core_Integrator_h

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
  class ParametersList;
  class Event;
  class GridParameters;
  namespace utils { class Timer; }
  /// Flavour of integration algorithm
  enum class IntegratorType {
    plain = 0, ///< Simple trial-and-error algorithm
    Vegas = 1, ///< Vegas algorithm (G.P. Lepage, 1977 \cite Lepage:1977sw)
    MISER = 2  ///< MISER stratified sampling algorithm
  };
  /// Monte-Carlo integrator instance
  class Integrator
  {
    public:
      /// Book the memory slots and structures for the integrator
      Integrator( const ParametersList& params );
      /// Class destructor
      ~Integrator();
      /**
       * Specify the function to be integrated
       * \param[in] ndim Number of dimensions on which the function will be integrated
       * \param[in] integrand Function to be integrated
       * \param[inout] params Run parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
       */
      void setFunction( unsigned int ndim, double integrand( double*, size_t, void* ), Parameters& params );
      /**
       * Algorithm to perform the n-dimensional Monte Carlo integration of a given function.
       * \param[out] result_ The cross section as integrated for the given phase space restrictions
       * \param[out] abserr_ The uncertainty associated to the computed cross section
       */
      virtual void integrate( double& result_, double& abserr_ ) = 0;
      /// Dimensional size of the phase space
      unsigned short dimensions() const;
      /// Generate a single event
      /// \param[in] callback The callback function applied on every event generated
      void generateOne( Event::callback callback = nullptr );
      /// Launch the event generation for a given number of events
      /// \param[in] callback The callback function applied on every event generated
      void generate( unsigned long num_events = 0, Event::callback callback = nullptr );

    protected:
      /// Compute the function value at the given phase space point
      virtual double eval( const std::vector<double>& x );
      unsigned int ncvg_; ///< Number of function calls to be computed for each point
      unsigned long seed_; ///< Random number generator seed
      struct gsl_rng_deleter
      {
        inline void operator()( gsl_rng* rng ) { gsl_rng_free( rng ); }
      };
      /// Instance of random number generator service
      std::unique_ptr<gsl_rng,gsl_rng_deleter> rng_;
      /// Set of parameters for the integration/event generation grid
      std::unique_ptr<GridParameters> grid_;
      IntegratorType type_; ///< Integration algorithm
      /// List of parameters to specify the integration range and the
      /// physics determining the phase space
      Parameters* input_params_;
      /// GSL structure storing the function to be integrated by this
      /// integrator instance (along with its parameters)
      std::unique_ptr<gsl_monte_function> function_;
      double result_, err_result_;

    private:
      /**
       * Store the event characterized by its _ndim-dimensional point in
       * the phase space to the output file
       * \brief Store the event in the output file
       * \param[in] x The d-dimensional point in the phase space defining the unique event to store
       * \param[in] callback The callback function for every event generated
       * \return A boolean stating whether or not the event could be saved
       */
      bool storeEvent( const std::vector<double>& x, Event::callback callback = nullptr );
      /// Start the correction cycle on the grid
      /// \param x Point in the phase space considered
      /// \param has_correction Correction cycle started?
      bool correctionCycle( std::vector<double>& x, bool& has_correction );
      /**
       * Set all the generation mode variables and align them to the
       *  integration grid set while computing the cross-section
       * \brief Prepare the class for events generation
       */
      void computeGenerationParameters();
      /// Generate a uniformly distributed (between 0 and 1) random number
      double uniform() const;
      /// Selected bin at which the function will be evaluated
      int ps_bin_;
      static constexpr int INVALID_BIN = -999;
  };
  std::ostream& operator<<( std::ostream&, const IntegratorType& );
}

#endif
