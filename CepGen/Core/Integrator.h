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
  struct GridParameters {
    GridParameters();
    unsigned int max;
    /// Has the generation been prepared using @a SetGen call? (very time-consuming operation, thus needs to be called once)
    bool gen_prepared;
    /// Maximal value of the function at one given point
    std::vector<double> f_max;
    /// Maximal value of the function in the considered integration range
    double f_max_global;
    std::vector<int> n;

    /// Maximal number of dimensions handled by this integrator instance
    static constexpr unsigned short max_dimensions_ = 15;
    /// Integration grid size parameter
    static constexpr unsigned short mbin_ = 3;
    static constexpr double inv_mbin_ = 1./mbin_;
  };
  /**
   * Main occurence of the Monte-Carlo integrator @cite PeterLepage1978192 developed by G.P. Lepage in 1978
   * \brief Monte-Carlo integrator instance
   */
  class Integrator
  {
    public:
      enum class Type { plain = 0, Vegas = 1, MISER = 2 };
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

      GridParameters grid;

    private:
      void computeGenerationParameters();
      double eval( const std::vector<double>& x );
      /// List of parameters to specify the integration range and the physics determining the phase space
      Parameters* input_params_;
      /// GSL structure storing the function to be integrated by this integrator instance (along with its parameters)
      std::unique_ptr<gsl_monte_function> function_;
      std::unique_ptr<gsl_rng,void(*)( gsl_rng* )> rng_;
      std::mutex mutex_;
  };
  std::ostream& operator<<( std::ostream&, const Integrator::Type& );
}

#endif

