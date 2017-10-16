#ifndef CepGen_Core_Integrator_h
#define CepGen_Core_Integrator_h

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_rng.h>

#include <vector>
#include <memory>

namespace CepGen
{
  class Parameters;
  /**
   * Main occurence of the Monte-Carlo integrator @cite PeterLepage1978192 developed by G.P. Lepage in 1978
   * \brief Monte-Carlo integrator instance
   */
  class Integrator
  {
    public:
      enum Type { Vegas = 1, MISER = 2 };
      /**
       * Book the memory slots and structures for the integrator
       * \note This code is based on the Vegas Monte Carlo integration algorithm developed by P. Lepage, as documented in @cite PeterLepage1978192
       * \param[in] dim_ Number of dimensions on which the function will be integrated
       * \param[in] f_ Function to be integrated
       * \param[inout] inParam_ Run parameters to define the phase space on which this integration is performed (embedded in an Parameters object)
       */
      Integrator( const unsigned int dim_, double f_(double*,size_t,void*), Parameters* inParam_, const Type& type = Vegas );
      /// Class destructor
      ~Integrator();
      /**
       * Algorithm to perform the n-dimensional Monte Carlo integration of a given function as described in @cite PeterLepage1978192
       * \author Primary author: G.P. Lepage
       * \author This C++ implementation: GSL
       * \param[out] result_ The cross section as integrated for the given phase space restrictions
       * \param[out] abserr_ The error associated to the computed cross section
       * \return 0 if the integration was performed successfully
       */
      int integrate( double& result_,double& abserr_ );
      /// Launch the generation of events
      void generate();
      /**
       * Generate one event according to the grid parameters set in \a SetGen
       * \brief Generate one single event according to the method defined in the Fortran 77 version of LPAIR
       * \return A boolean stating if the generation was successful (in term of the computed weight for the phase space point)
       */
      bool generateOneEvent();
      /// Dimensional size of the phase space
      const unsigned short dimensions() const { return ( !function_ ) ? 0 : function_->dim; }
    private:
      /**
       * Evaluate the function to be integrated at a point @a x_, using the default Parameters object @a fInputParameters
       * \param[in] x_ The point at which the function is to be evaluated
       * \return Function value at this point @a x_
       */
      inline double F( const std::vector<double>& x ) { return F( x, input_params_ ); }
      /**
       * Evaluate the function to be integrated at a point @a x_, given a set of Parameters @a ip_
       * \param[in] x_ The point at which the function is to be evaluated
       * \param[in] ip_ A set of parameters to fully define the function
       * \return Function value at this point \a x
       */
      inline double F( const std::vector<double>& x, Parameters* ip ) {
        return function_->f( (double*)&x[0], function_->dim, (void*)ip );
      }
      /**
       * Store the event characterized by its _ndim-dimensional point in the phase
       * space to the output file
       * \brief Store the event in the output file
       * \param[in] x The d-dimensional point in the phase space defining the unique event to store
       * \return A boolean stating whether or not the event could be saved
       */
      bool storeEvent( const std::vector<double>& x );
      /// Start the correction cycle on the grid
      /// \param x Point in the phase space considered
      /// \param has_correction Correction cycle started?
      bool correctionCycle( std::vector<double>& x, bool& has_correction );
      /**
       * Set all the generation mode variables and align them to the integration grid set while computing the cross-section
       * \brief Prepare the class for events generation
       */
      void setGen();
      double uniform() const { return gsl_rng_uniform( rng_ ); }

      /// Maximal number of dimensions handled by this integrator instance
      static constexpr unsigned short max_dimensions_ = 15;
      /// Integration grid size parameter
      static constexpr unsigned short mbin_ = 3;
      static constexpr double inv_mbin_ = 1./mbin_;

      Type algorithm_;
      /// Selected bin at which the function will be evaluated
      int ps_bin_;
      double correc_;
      double correc2_;
      /// List of parameters to specify the integration range and the physics determining the phase space
      Parameters* input_params_;
      /// Has the grid been prepared for integration?
      bool grid_prepared_;
      /// Has the generation been prepared using @a SetGen call? (very time-consuming operation, thus needs to be called once)
      bool gen_prepared_;
      /// Maximal value of the function at one given point
      std::vector<double> f_max_;
      double f_max2_;
      double f_max_diff_;
      double f_max_old_;
      /// Maximal value of the function in the considered integration range
      double f_max_global_;
      std::vector<int> n_;
      std::vector<int> nm_;
      /// GSL structure storing the function to be integrated by this integrator instance (along with its parameters)
      std::unique_ptr<gsl_monte_function> function_;
      gsl_rng* rng_;
      /// Number of function calls to be computed for each point
      int num_converg_;
      /// Number of iterations for the integration
      unsigned int num_iter_;
  };
}

#endif

