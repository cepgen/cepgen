#ifndef CepGen_Core_GeneratorWorker_h
#define CepGen_Core_GeneratorWorker_h

#include "CepGen/Event/Event.h"

#include <vector>
#include <memory>

namespace cepgen
{
  class Integrator;
  class Integrand;
  class Parameters;
  class GridParameters;
  /// Monte-Carlo generator instance
  class GeneratorWorker
  {
    public:
      /// Book the memory slots and structures for the generator
      explicit GeneratorWorker( Parameters* );
      /// Specify the integrator instance handled by the mother generator
      void setIntegrator( const Integrator* integ );
      /// Launch the event generation
      /// \param[in] num_events Number of events to generate
      /// \param[in] callback The callback function applied on every event generated
      void generate( size_t num_events = 0, Event::callback callback = nullptr );
      /// Function evaluator
      Integrand& integrand() { return *integrand_; }

    private:
      /// Generate a single event
      /// \param[in] callback The callback function applied on every event generated
      bool next( Event::callback callback = nullptr );
      /// Store the event in the output file
      /// \param[in] callback The callback function for every event generated
      /// \return A boolean stating whether or not the event was successfully saved
      bool storeEvent( Event::callback );
      /// Apply a correction cycle to the grid
      bool correctionCycle( bool& );
      /// Prepare the object for event generation
      void computeGenerationParameters();

      /// Local event weight evaluator
      std::unique_ptr<Integrand> integrand_;
      /// Pointer to the mother-handled integrator instance
      const Integrator* integrator_; // not owning
      /// Steering parameters for the event generation
      const Parameters* params_; // not owning
      /// Set of parameters for the integration/event generation grid
      std::unique_ptr<GridParameters> grid_;
      /// Selected bin at which the function will be evaluated
      int ps_bin_; ///< Last bin to be corrected
      std::vector<double> coords_; ///< Phase space coordinates being evaluated
      static constexpr int INVALID_BIN = -999; ///< Placeholder for invalid bin indexing
  };
}

#endif
