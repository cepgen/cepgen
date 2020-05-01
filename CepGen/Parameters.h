#ifndef CepGen_Parameters_h
#define CepGen_Parameters_h

#include "CepGen/Physics/Kinematics.h"

#include <memory>

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>

namespace cepgen
{
  class EventModifier;
  class ParametersList;
  namespace proc { class Process; }
  namespace io { class ExportModule; }
  namespace utils {
    class TamingFunction;
    class TimeKeeper;
  }
  enum class IntegratorType;
  /// An ordered collection of event modification algorithms
  typedef std::vector<std::unique_ptr<EventModifier> > EventModifiersSequence;
  /// An ordered collection of event export modules
  typedef std::vector<std::unique_ptr<io::ExportModule> > ExportModulesSequence;
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
      /// Dump the input parameters in the terminal
      friend std::ostream& operator<<( std::ostream&, const Parameters* );

      /// Initialise the timekeeper instance
      void setTimeKeeper( utils::TimeKeeper* );
      /// Pointer to a timekeeper instance
      utils::TimeKeeper* timeKeeeper() { return tmr_.get(); }

      /// Common user-defined parameters
      std::shared_ptr<ParametersList> general;
      /// Integrator specific user-defined parameters
      std::shared_ptr<ParametersList> integrator;

      //----- process to compute

      /// Is this parameters collection holding any physics process?
      bool hasProcess() const { return !( !process_ ); }
      /// Process for which the cross-section will be computed and the events will be generated
      proc::Process& process();
      /// Process for which the cross-section will be computed and the events will be generated
      const proc::Process& process() const;
      /// Name of the process considered
      std::string processName() const;
      /// Remove the process pointer
      void clearProcess();
      /// Copy a process configuration
      void setProcess( std::unique_ptr<proc::Process> proc );
      /// Set a process configuration
      void setProcess( proc::Process* proc );

      //----- events kinematics

      /// Events kinematics for phase space definition
      Kinematics kinematics;

      //----- events generation

      /// Collection of events generation parameters
      struct Generation
      {
        Generation();
        Generation( const Generation& );
        bool enabled; ///< Are we generating events ? (true) or are we only computing the cross-section ? (false)
        unsigned long maxgen; ///< Maximal number of events to generate in this run
        bool symmetrise; ///< Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
        unsigned int gen_print_every; ///< Frequency at which the events are displayed to the end-user
        unsigned int num_threads; ///< Number of threads to perform the integration
        unsigned int num_points; ///< Number of points to "shoot" in each integration bin by the algorithm
      };
      /// Get the events generation parameters
      Generation& generation() { return generation_; }
      /// Get the events generation parameters
      const Generation& generation() const { return generation_; }

      //----- event modification (e.g. hadronisation, decay) algorithm

      /// Event modification algorithm to use
      EventModifier& eventModifier( size_t );
      /// Retrieve the list of event modification algorithms to run
      EventModifiersSequence& eventModifiersSequence() { return evt_modifiers_; }
      /// Retrieve the list of event modification algorithms to run
      const EventModifiersSequence& eventModifiersSequence() const { return evt_modifiers_; }
      /// Add a new event modification algorithm to the sequence
      void addModifier( std::unique_ptr<EventModifier> );
      /// Add a new event modification algorithm to the sequence
      void addModifier( EventModifier* );

      //----- event output algorithms

      /// Output module
      io::ExportModule& outputModule( size_t );
      /// Retrieve the list of output modules to run
      ExportModulesSequence& outputModulesSequence() { return out_modules_; }
      /// Retrieve the list of output modules to run
      const ExportModulesSequence& outputModulesSequence() const { return out_modules_; }
      /// Set a new output module definition
      void addOutputModule( std::unique_ptr<io::ExportModule> mod );
      /// Set the pointer to a output module
      void addOutputModule( io::ExportModule* mod );

      //----- taming functions

      /// Functionals to be used to account for rescattering corrections (implemented within the process)
      std::vector<utils::TamingFunction> taming_functions;

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
      /// Physics process held by these parameters
      std::unique_ptr<proc::Process> process_;
      /// Collection of event modification algorithms to be applied
      EventModifiersSequence evt_modifiers_;
      /// Collection of event output modules to be applied
      ExportModulesSequence out_modules_;
      /// Total generation time (in seconds)
      double total_gen_time_;
      /// Number of events already generated
      unsigned long num_gen_events_;
      /// Events generation parameters
      Generation generation_;
      /// A collection of stopwatches for timing
      std::unique_ptr<utils::TimeKeeper> tmr_;
  };
}

#endif
