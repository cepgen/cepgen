#ifndef CepGen_Parameters_h
#define CepGen_Parameters_h

#include "CepGen/Physics/Kinematics.h"

#include <memory>

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>

namespace cepgen
{
  class Event;
  class EventModifier;
  class ParametersList;
  namespace proc { class GenericProcess; }
  namespace io { class GenericExportHandler; }
  namespace utils { class TamingFunctionsCollection; }
  enum class IntegratorType;
  typedef std::vector<std::unique_ptr<EventModifier> > EventModifiersSequence;
  /// List of parameters used to start and run the simulation job
  class Parameters
  {
    public:
      Parameters();
      /// Copy constructor (transfers ownership to the process/event modification algorithm!)
      Parameters( Parameters& );
      /// Const copy constructor (all but the process and the event modification algorithm)
      Parameters( const Parameters& );
      ~Parameters(); // required for unique_ptr initialisation!

      /// Assignment operator
      Parameters& operator=( Parameters );

      /// Set the polar angle range for the produced leptons
      /// \param[in] thetamin The minimal value of \f$\theta\f$ for the outgoing leptons
      /// \param[in] thetamax The maximal value of \f$\theta\f$ for the outgoing leptons
      void setThetaRange( float thetamin, float thetamax );
      /// Dump the input parameters in the terminal
      friend std::ostream& operator<<( std::ostream&, const Parameters* );

      std::shared_ptr<ParametersList> general;

      //----- process to compute

      /// Process for which the cross-section will be computed and the events will be generated
      proc::GenericProcess* process();
      /// Process for which the cross-section will be computed and the events will be generated
      const proc::GenericProcess* process() const;
      /// Name of the process considered
      std::string processName() const;
      /// Set the process to study
      void setProcess( std::unique_ptr<proc::GenericProcess> proc );
      /// Set the process to study
      void setProcess( proc::GenericProcess* proc );

      //----- events kinematics

      /// Events kinematics for phase space definition
      Kinematics kinematics;

      //----- VEGAS

      /// Collection of integrator parameters
      struct Integration
      {
        Integration();
        Integration( const Integration& );
        ~Integration();
        IntegratorType type;
        unsigned int ncvg; ///< Number of function calls to be computed for each point
        long rng_seed; ///< Random number generator seed
        gsl_rng_type* rng_engine; ///< Random number generator engine
        gsl_monte_vegas_params vegas;
        double vegas_chisq_cut;
        gsl_monte_miser_params miser;
        double result, err_result;
      };
      Integration& integration() { return integration_; }
      const Integration& integration() const { return integration_; }

      //----- events generation

      /// Collection of events generation parameters
      struct Generation
      {
        Generation();
        Generation( const Generation& );
        bool enabled; ///< Are we generating events ? (true) or are we only computing the cross-section ? (false)
        unsigned long maxgen; ///< Maximal number of events to generate in this run
        bool symmetrise; ///< Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
        bool treat; ///< Is the integrand to be smoothed for events generation?
        unsigned int gen_print_every; ///< Frequency at which the events are displayed to the end-user
        unsigned int num_threads; ///< Number of threads to perform the integration
        unsigned int num_points; ///< Number of points to "shoot" in each integration bin by the algorithm
      };
      Generation& generation() { return generation_; }
      const Generation& generation() const { return generation_; }

      /// Specify if the generated events are to be stored
      void setStorage( bool store ) { store_ = store; }
      /// Are the events generated in this run to be stored in the output file ?
      bool storage() const { return store_; }

      /// Set a new output module definition
      void setOutputModule( std::unique_ptr<io::GenericExportHandler> mod );
      /// Set the pointer to a output module
      void setOutputModule( io::GenericExportHandler* mod );
      /// Output module definition
      io::GenericExportHandler* outputModule();

      //----- event modification (e.g. hadronisation, decay) algorithm

      /// Event modification algorithm to use
      EventModifier* eventModifier( size_t );
      /// Retrieve the list of event modification algorithms to run
      EventModifiersSequence& eventModifiersSequence() { return evt_modifiers_; }
      /// Retrieve the list of event modification algorithms to run
      const EventModifiersSequence& eventModifiersSequence() const { return evt_modifiers_; }
      /// Name of the modification algorithm (if applicable)
      std::string eventModifierName( size_t ) const;
      /// Add a new event modification algorithm to the sequence
      void addModifier( std::unique_ptr<EventModifier> );
      /// Add a new event modification algorithm to the sequence
      void addModifier( EventModifier* );
      /// Set the event modification algorithms sequence
      void setModifiersSequence( EventModifiersSequence& );

      //----- taming functions

      /// Functionals to be used to account for rescattering corrections (implemented within the process)
      std::shared_ptr<utils::TamingFunctionsCollection> taming_functions;

      //----- run operations

      /// Reset the total generation time and the number of events generated for this run, prepare kinematics
      void prepareRun();
      /// Add a new timing into the total generation time
      /// \param[in] gen_time Time to add (in seconds)
      void addGenerationTime( double gen_time );
      /// Return the total generation time for this run (in seconds)
      inline double totalGenerationTime() const { return total_gen_time_; }
      /// Total number of events already generated in this run
      inline unsigned int numGeneratedEvents() const { return num_gen_events_; }

    private:
      std::unique_ptr<proc::GenericProcess> process_;
      EventModifiersSequence evt_modifiers_;
      /// Storage object
      std::unique_ptr<io::GenericExportHandler> out_module_;

      bool store_;
      /// Total generation time (in seconds)
      double total_gen_time_;
      /// Number of events already generated
      unsigned long num_gen_events_;
      /// Integrator parameters
      Integration integration_;
      /// Events generation parameters
      Generation generation_;
  };
}

#endif
