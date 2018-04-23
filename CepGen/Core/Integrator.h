#ifndef CepGen_Core_Integrator_h
#define CepGen_Core_Integrator_h

#include <vector>
#include <memory>
#include <mutex>
#include <functional>

#include <string.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

namespace CepGen
{
  class Parameters;
  class Event;
  class GridParameters;
  /**
   * Main occurence of the Monte-Carlo integrator @cite PeterLepage1978192 developed by G.P. Lepage in 1978
   * \brief Monte-Carlo integrator instance
   */
  class Integrator
  {
    public:
      enum class Type { plain = 0, Vegas = 1, MISER = 2 };
      enum class VegasMode { importance = 1, importanceOnly = 0, stratified = -1 };
      /**
       * Book the memory slots and structures for the integrator
       * \note Three integration algorithms are currently supported:
       *  * the plain algorithm randomly sampling points in the phase space
       *  * the Vegas algorithm developed by P. Lepage, as documented in @cite PeterLepage1978192,
       *  * the MISER algorithm developed by W.H. Press and G.R. Farrar, as documented in @cite Press:1989vk.
       * \param[in] ndim Number of dimensions on which the function will be integrated
       * \param[in] integrand Function to be integrated
       * \param[inout] params Run parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
       */
      Integrator( unsigned int ndim, double integrand(double*,size_t,void*), Parameters* params );
      /// Class destructor
      ~Integrator();
      /**
       * Algorithm to perform the n-dimensional Monte Carlo integration of a given function.
       * \author This C++ implementation: GSL
       * \param[out] result_ The cross section as integrated for the given phase space restrictions
       * \param[out] abserr_ The error associated to the computed cross section
       * \return 0 if the integration was performed successfully
       */
      int integrate( double& result_, double& abserr_ );
      /// Dimensional size of the phase space
      unsigned short dimensions() const;
      void generateOne( std::function<void( const Event&, unsigned long )> callback = nullptr );
      void generate( unsigned long num_events = 0, std::function<void( const Event&, unsigned long )> callback = nullptr );

    private:
      void computeGenerationParameters();
      double uniform() const;
      double eval( const std::vector<double>& x );
      /// List of parameters to specify the integration range and the physics determining the phase space
      Parameters* input_params_;
      /// GSL structure storing the function to be integrated by this integrator instance (along with its parameters)
      std::unique_ptr<gsl_monte_function> function_;
      std::unique_ptr<gsl_rng,void(*)( gsl_rng* )> rng_;
      std::unique_ptr<GridParameters> grid_;
      std::mutex mutex_;
  };
  std::ostream& operator<<( std::ostream&, const Integrator::Type& );
  std::ostream& operator<<( std::ostream&, const Integrator::VegasMode& );
}

#endif

