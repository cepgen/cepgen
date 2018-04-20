#ifndef CepGen_Parameters_h
#define CepGen_Parameters_h

#include "CepGen/Core/Integrator.h"
#include "CepGen/Physics/Kinematics.h"

#include <memory>

#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>

namespace CepGen
{
  class Event;
  class TamingFunctionsCollection;
  namespace Process { class GenericProcess; }
  namespace Hadroniser { class GenericHadroniser; }
  /// List of parameters used to start and run the simulation job
  class Parameters
  {
    public:
      Parameters();
      /// Copy constructor (transfers ownership to the process/hadroniser!)
      Parameters( Parameters& );
      /// Const copy constructor (all but the process and the hadroniser)
      Parameters( const Parameters& );
      ~Parameters(); // required for unique_ptr initialisation!
      /// Set the polar angle range for the produced leptons
      /// \param[in] thetamin The minimal value of \f$\theta\f$ for the outgoing leptons
      /// \param[in] thetamax The maximal value of \f$\theta\f$ for the outgoing leptons
      void setThetaRange( float thetamin, float thetamax );
      /// Dump the input parameters in the console
      void dump( std::ostream& os = *Logger::get().output, bool pretty = true ) const;

      //----- process to compute

      /// Process for which the cross-section will be computed and the events will be generated
      Process::GenericProcess* process();
      /// Name of the process considered
      std::string processName() const;
      /// Set the process to study
      void setProcess( Process::GenericProcess* proc );
      std::unique_ptr<Process::GenericProcess> processClone() const;

      //----- events kinematics

      /// Events kinematics for phase space definition
      Kinematics kinematics;

      //----- VEGAS

      /// Collection of integrator parameters
      struct IntegratorParameters
      {
        IntegratorParameters();
        Integrator::Type type;
        /// Number of function calls to be computed for each point
        unsigned int ncvg; // ??
        /// Number of points to "shoot" in each integration bin by the algorithm
        unsigned int npoints;
        /// Random number generator seed
        long rng_seed;
        /// Random number generator engine
        gsl_rng_type* rng_engine;
        gsl_monte_vegas_params vegas;
        double vegas_chisq_cut;
        gsl_monte_miser_params miser;
        double result, err_result;
      };
      /// Integrator parameters
      IntegratorParameters integrator;

      //----- events generation

      /// Collection of events generation parameters
      struct Generation
      {
        Generation();
        /// Are we generating events ? (true) or are we only computing the cross-section ? (false)
        bool enabled;
        /// Maximal number of events to generate in this run
        unsigned int maxgen;
        /// Do we want the events to be symmetrised with respect to the \f$z\f$-axis ?
        bool symmetrise;
        /// Number of events already generated in this run
        unsigned int ngen;
        /// Frequency at which the events are displayed to the end-user
        unsigned int gen_print_every;
        /// Number of threads to perform the integration
        unsigned int num_threads;
      };
      /// Events generation parameters
      Generation generation;

      /// Specify if the generated events are to be stored
      void setStorage( bool store ) { store_ = store; }
      /// Are the events generated in this run to be stored in the output file ?
      bool storage() const { return store_; }

      //----- hadronisation algorithm

      /// Hadronisation algorithm to use for the proton(s) fragmentation
      Hadroniser::GenericHadroniser* hadroniser();
      /// Set the hadronisation algorithm
      void setHadroniser( Hadroniser::GenericHadroniser* hadr );
      /// Maximal number of trials for the hadronisation of the proton(s) remnants
      unsigned int hadroniser_max_trials;

      //----- taming functions

      /// Functionals to be used to account for rescattering corrections (implemented within the process)
      std::shared_ptr<TamingFunctionsCollection> taming_functions;

      //----- run statistics

      /// Reset the total generation time and the number of events generated for this run
      void clearRunStatistics();
      /// Add a new timing into the total generation time
      /// \param[in] gen_time Time to add (in seconds)
      void addGenerationTime( double gen_time );
      /// Return the total generation time for this run (in seconds)
      inline double totalGenerationTime() const { return total_gen_time_; }
      /// Total number of events already generated in this run
      inline unsigned int numGeneratedEvents() const { return num_gen_events_; }

    private:
      std::unique_ptr<Process::GenericProcess> process_;
      std::shared_ptr<Hadroniser::GenericHadroniser> hadroniser_;

      bool store_;
      /// Total generation time (in seconds)
      double total_gen_time_;
      /// Number of events already generated
      unsigned int num_gen_events_;
  };
}

#endif
