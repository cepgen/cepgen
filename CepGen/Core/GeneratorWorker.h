#ifndef CepGen_Core_GeneratorWorker_h
#define CepGen_Core_GeneratorWorker_h

#include "CepGen/Event/Event.h"

#include <vector>
#include <memory>

namespace cepgen
{
  class Integrator;
  class Parameters;
  class GridParameters;
  /// Monte-Carlo generator instance
  class GeneratorWorker
  {
    public:
      /// Book the memory slots and structures for the generator
      explicit GeneratorWorker( Parameters* );
      /// Specify the integrator instance handled by the mother generator
      void setIntegrator( Integrator* integ );
      /// Launch the event generation
      /// \param[in] num_events Number of events to generate
      /// \param[in] callback The callback function applied on every event generated
      void generate( size_t num_events = 0, Event::callback callback = nullptr );

    private:
      /// Generate a single event
      /// \param[in] callback The callback function applied on every event generated
      void generateOne( Event::callback callback = nullptr );
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

      /// Steering parameters for the event generation
      Parameters* input_params_; // not owning
      /// Pointer to the mother-handled integrator instance
      Integrator* integrator_; // not owning
      /// Set of parameters for the integration/event generation grid
      std::unique_ptr<GridParameters> grid_;
      /// Selected bin at which the function will be evaluated
      int ps_bin_;
      bool initialised_ = false;
      static constexpr int INVALID_BIN = -999;
  };
}

#endif
