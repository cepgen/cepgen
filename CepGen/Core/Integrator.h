#ifndef CepGen_Core_Integrator_h
#define CepGen_Core_Integrator_h

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_rng.h>

#include <vector>
#include <memory>
#include <functional>

namespace CepGen
{
  class Parameters;
  class Event;
  struct GridParameters {
    GridParameters() :
      grid_prepared( false ), gen_prepared( false ),
      correc( 0. ), correc2( 0. ),
      f_max2( 0. ), f_max_diff( 0. ), f_max_old( 0. ), f_max_global( 0. ) {}
    /// Has the grid been prepared for integration?
    bool grid_prepared;
    /// Has the generation been prepared using @a SetGen call? (very time-consuming operation, thus needs to be called once)
    bool gen_prepared;
    double correc;
    double correc2;
    /// Maximal value of the function at one given point
    std::vector<double> f_max;
    double f_max2;
    double f_max_diff;
    double f_max_old;
    /// Maximal value of the function in the considered integration range
    double f_max_global;
    std::vector<int> n;
    std::vector<int> nm;

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
      enum Type { Plain = 0, Vegas = 1, MISER = 2 };
      /**
       * Book the memory slots and structures for the integrator
       * \note Three integration algorithms are currently supported:
       *  * the plain algorithm randomly sampling points in the phase space
       *  * the Vegas algorithm developed by P. Lepage, as documented in @cite PeterLepage1978192,
       *  * the MISER algorithm developed by W.H. Press and G.R. Farrar, as documented in @cite Press:1989vk.
       * \param[in] dim_ Number of dimensions on which the function will be integrated
       * \param[in] f_ Function to be integrated
       * \param[inout] inParam_ Run parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
       */
      Integrator( const unsigned int dim_, double f_(double*,size_t,void*), Parameters* inParam_ );
      /// Class destructor
      ~Integrator();
      /**
       * Algorithm to perform the n-dimensional Monte Carlo integration of a given function.
       * \author This C++ implementation: GSL
       * \param[out] result_ The cross section as integrated for the given phase space restrictions
       * \param[out] abserr_ The error associated to the computed cross section
       * \return 0 if the integration was performed successfully
       */
      int integrate( double& result_,double& abserr_ );
      /// Dimensional size of the phase space
      unsigned short dimensions() const;
      void generate( unsigned long num_events, std::function<void( const Event&, unsigned long )> callback = nullptr );

    private:
      void setGen();
      /// List of parameters to specify the integration range and the physics determining the phase space
      Parameters* input_params_;
      GridParameters grid_;
      /// GSL structure storing the function to be integrated by this integrator instance (along with its parameters)
      std::unique_ptr<gsl_monte_function> function_;
      std::shared_ptr<gsl_rng> rng_;
  };
  std::ostream& operator<<( std::ostream&, const Integrator::Type& );

  class ThreadWorker
  {
    public:
      ThreadWorker( std::shared_ptr<gsl_rng> rng, gsl_monte_function* function, GridParameters* grid, std::function<void( const Event&, unsigned long )> callback = nullptr );

      /// Generate one event according to the grid parameters set in the initialisation
      /// \return A boolean stating if the generation was successful (in term of the computed weight for the phase space point)
      bool generate();

    private:
      bool next();
      double uniform() const;
      /// Start the correction cycle on the grid
      /// \param x Point in the phase space considered
      /// \param has_correction Correction cycle started?
      bool correctionCycle( std::vector<double>& x, bool& has_correction );
      /// Evaluate the function to be integrated at a point
      /// \param[in] x The point at which the function is to be evaluated
      /// \return Function value at this point
      double eval( const std::vector<double>& x );
      /// Store the event in the output file
      /// \param[in] x The d-dimensional point in the phase space defining the unique event to store
      /// \return A boolean stating whether or not the event could be saved
      bool storeEvent( const std::vector<double>& x );

      /// Selected bin at which the function will be evaluated
      int ps_bin_;

      std::shared_ptr<gsl_rng> rng_;
      std::shared_ptr<gsl_monte_function> function_;
      std::shared_ptr<GridParameters> grid_;

      Parameters* params_;
      std::function<void( const Event&, unsigned long )> callback_;
  };
}

#endif

