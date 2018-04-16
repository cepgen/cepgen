#ifndef CepGen_Core_ThreadWorker_h
#define CepGen_Core_ThreadWorker_h

#include <vector>
#include <memory>
#include <mutex>
#include <functional>

#include "CepGen/Processes/GenericProcess.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_rng.h>

namespace CepGen
{
  //--- forward declarations
  class Parameters;
  class Event;
  class GridParameters;

  //--- class definition
  /// Single thread handling events definition for a given parameters set
  class ThreadWorker
  {
    public:
      /// Build a thread using a Mutex, GSL objects, a grid and the functional
      ThreadWorker( std::mutex* mutex, gsl_rng* rng, gsl_monte_function* function, GridParameters* grid, std::function<void( const Event&, unsigned long )>& callback );

      /// Generate one event according to the grid parameters set in the initialisation
      /// \return A boolean stating if the generation was successful (in term of the computed weight for the phase space point)
      bool generate();

    private:
      /// Generate the next event
      bool next();
      /// Generate a random number following a uniform distribution
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

      /// Pointer to the random numbers generator
      std::shared_ptr<gsl_rng> rng_;
      /// Pointer to the functional to be used for generations
      gsl_monte_function* function_;

      GridParameters* grid_;
      std::vector<int> grid_nm_;
      double grid_correc_;
      double grid_correc2_;
      double grid_f_max2_;
      double grid_f_max_diff_;
      double grid_f_max_old_;

      /// Pointer to the (thread-independent) run parameters
      Parameters* global_params_;
      /// Pointer to the run parameters specific to this thread
      Parameters* local_params_;
      /// Clone of the user physics process to be handled by this thread
      std::unique_ptr<Process::GenericProcess> process_;
      std::mutex* mutex_;

      /// A function to be called after each event generation
      std::function<void( const Event&, unsigned long )> callback_;
  };
}

#endif
