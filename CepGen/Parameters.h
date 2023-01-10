/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CepGen_Parameters_h
#define CepGen_Parameters_h

#include <memory>

#include "CepGen/Physics/Kinematics.h"

namespace cepgen {
  class EventExporter;
  class EventModifier;
  class ParametersList;
  namespace proc {
    class Process;
  }
  namespace utils {
    class TimeKeeper;
    class Functional;
  }  // namespace utils
  enum class IntegratorType;
  /// An ordered collection of event modification algorithms
  typedef std::vector<std::unique_ptr<EventModifier> > EventModifiersSequence;
  /// An ordered collection of event export modules
  typedef std::vector<std::unique_ptr<EventExporter> > EventExportersSequence;
  /// An ordered collection of taming functions evaluators
  typedef std::vector<std::unique_ptr<utils::Functional> > TamingFunctionsSequence;
  /// List of parameters used to start and run the simulation job
  class Parameters {
  public:
    Parameters();
    /// Copy constructor (transfers ownership to the process/event modification algorithm!)
    Parameters(Parameters&);
    /// Const copy constructor (all but the process and the event modification algorithm)
    Parameters(const Parameters&);
    ~Parameters();  // required for unique_ptr initialisation!

    /// Assignment operator
    Parameters& operator=(Parameters);
    /// Dump the input parameters in the terminal
    friend std::ostream& operator<<(std::ostream&, const Parameters*);

    /// Initialise the timekeeper instance
    void setTimeKeeper(utils::TimeKeeper*);
    /// Pointer to a timekeeper instance
    utils::TimeKeeper* timeKeeper() { return tmr_.get(); }
    /// Pointer to a timekeeper instance
    const utils::TimeKeeper* timeKeeper() const { return tmr_.get(); }

    /// Integrator specific user-defined parameters
    ParametersList par_integrator;

    //----- process to compute

    /// Is this parameters collection holding any physics process?
    bool hasProcess() const { return !(!process_); }
    /// Process for which the cross-section will be computed and the events will be generated
    proc::Process& process();
    /// Process for which the cross-section will be computed and the events will be generated
    const proc::Process& process() const;
    /// Name of the process considered
    std::string processName() const;
    /// Remove the process pointer
    void clearProcess();
    /// Copy a process configuration
    void setProcess(std::unique_ptr<proc::Process> proc);
    /// Set a process configuration
    void setProcess(proc::Process* proc);

    //----- events kinematics

    /// Events kinematics for phase space definition
    const Kinematics& kinematics() const;

    //----- events generation

    /// Collection of events generation parameters
    class Generation : public SteeredObject<Generation> {
    public:
      /// Build a generation parameters collection from a user input
      explicit Generation(const ParametersList& = ParametersList());

      static ParametersDescription description();

      /// Set the target luminosity to reach (in pb^-1)
      void setTargetLuminosity(double lumi_invpb) { target_lumi_ = lumi_invpb; }
      /// Target luminosity to reach (in pb^-1)
      double targetLuminosity() const { return target_lumi_; }
      /// Set the maximal number of events to generate
      void setMaxGen(size_t max_gen) { max_gen_ = max_gen; }
      /// Maximal number of events to generate
      size_t maxGen() const { return max_gen_; }
      /// Are we generating events? (true) or only computing the cross-section? (false)
      bool enabled() const { return max_gen_ > 0; }
      /// Set the frequency at which events are displayed to the end-user
      void setPrintEvery(size_t print_every) { gen_print_every_ = print_every; }
      /// Frequency at which events are displayed to the end-user
      size_t printEvery() const { return gen_print_every_; }
      /// Switch on/off the symmetrisation of the z-axis for each event
      void setSymmetrise(bool sym) { symmetrise_ = sym; }
      /// Do we want the events to be symmetric with respect to the \f$z\f$-axis ?
      bool symmetrise() const { return symmetrise_; }
      /// Set the number of threads for the events generation
      void setNumThreads(size_t nt) { num_threads_ = nt; }
      /// Number of threads to perform the events generation
      size_t numThreads() const { return num_threads_; }
      /// Set the number of points to probe in each integration bin
      void setNumPoints(size_t np) { num_points_ = np; }
      /// Number of points to "shoot" in each integration bin by the algorithm
      size_t numPoints() const { return num_points_; }

    private:
      int max_gen_, gen_print_every_;
      double target_lumi_;
      bool symmetrise_;
      int num_threads_, num_points_;
    };
    /// Get the events generation parameters
    Generation& generation() { return generation_; }
    /// Get the events generation parameters
    const Generation& generation() const { return generation_; }

    //----- event modification (e.g. hadronisation, decay) algorithm

    /// Event modification algorithm to use
    EventModifier& eventModifier(size_t);
    /// Retrieve the list of event modification algorithms to run
    EventModifiersSequence& eventModifiersSequence() { return evt_modifiers_; }
    /// Retrieve the list of event modification algorithms to run
    const EventModifiersSequence& eventModifiersSequence() const { return evt_modifiers_; }
    /// Remove all event modifiers from sequence
    void clearEventModifiersSequence();
    /// Add a new event modification algorithm to the sequence
    void addModifier(std::unique_ptr<EventModifier>);
    /// Add a new event modification algorithm to the sequence
    void addModifier(EventModifier*);

    //----- event output algorithms

    /// Output module
    EventExporter& eventExporter(size_t);
    /// Retrieve the list of output modules to run
    EventExportersSequence& eventExportersSequence() { return evt_exporters_; }
    /// Retrieve the list of output modules to run
    const EventExportersSequence& eventExportersSequence() const { return evt_exporters_; }
    /// Remove all output modules from sequence
    void clearEventExportersSequence();
    /// Set a new output module definition
    void addEventExporter(std::unique_ptr<EventExporter> mod);
    /// Set the pointer to a output module
    void addEventExporter(EventExporter* mod);

    //----- taming functions

    /// List of all taming functions definitions
    const TamingFunctionsSequence& tamingFunctions() const { return taming_functions_; }
    /// Set a new taming function definition
    void addTamingFunction(std::unique_ptr<utils::Functional>);

    //----- run operations

    /// Reset the total generation time and the number of events generated for this run, prepare kinematics
    void prepareRun();
    /// Add a new timing into the total generation time
    /// \param[in] gen_time Time to add (in seconds)
    void addGenerationTime(double gen_time);
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
    EventExportersSequence evt_exporters_;
    /// Functions to be used to account for rescattering corrections
    TamingFunctionsSequence taming_functions_;
    /// Total generation time (in seconds)
    double total_gen_time_{0.};
    /// Number of events already generated
    unsigned long num_gen_events_{0ul};
    /// Events generation parameters
    Generation generation_;
    /// A collection of stopwatches for timing
    std::unique_ptr<utils::TimeKeeper> tmr_;
  };
}  // namespace cepgen

#endif
