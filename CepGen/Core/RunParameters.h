/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2025  Laurent Forthomme
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

#ifndef CepGen_Core_RunParameters_h
#define CepGen_Core_RunParameters_h

#include "CepGen/Physics/Kinematics.h"

namespace cepgen {
  class EventExporter;
  class EventModifier;
  class ParametersList;
  enum class IntegratorType;
}  // namespace cepgen
namespace cepgen::proc {
  class Process;
}  // namespace cepgen::proc
namespace cepgen::utils {
  class Functional;
  class TimeKeeper;
}  // namespace cepgen::utils

namespace cepgen {
  typedef std::vector<std::unique_ptr<EventModifier> > EventModifiersSequence;  ///< Event modification algos ordered set
  typedef std::vector<std::unique_ptr<EventExporter> > EventExportersSequence;  ///< Event export modules ordered set
  typedef std::vector<std::unique_ptr<utils::Functional> > TamingFunctionsSequence;  ///< Taming functions evaluators set
}  // namespace cepgen

namespace cepgen {
  /// List of parameters used to start and run the simulation job
  class RunParameters : public SteeredObject<RunParameters> {
  public:
    RunParameters();
    RunParameters(RunParameters&);  ///< Copy constructor (transfers ownership to process/event modification algorithm!)
    RunParameters(const RunParameters&);  ///< Const copy constructor (all but process + event handling algorithms)
    ~RunParameters() override;            // required for unique_ptr initialisation!

    static ParametersDescription description();

    RunParameters& operator=(RunParameters);  ///< Assignment operator

    friend std::ostream& operator<<(std::ostream&, const RunParameters&);  ///< User-readable dump of runtime parameters

    void setTimeKeeper(utils::TimeKeeper*);                             ///< Initialise the timekeeper instance
    utils::TimeKeeper* timeKeeper() { return timer_.get(); }              ///< Pointer to a timekeeper instance
    const utils::TimeKeeper* timeKeeper() const { return timer_.get(); }  ///< Pointer to a timekeeper instance

    void initialiseModules();  ///< Initialise the event handling modules for an event generation

    inline ParametersList& integrator() { return integrator_; }              ///< Integrator specific user parameters
    inline const ParametersList& integrator() const { return integrator_; }  ///< Integrator specific user parameters

    //----- process to compute

    inline bool hasProcess() const { return !(!process_); }  ///< Are we holding any physics process?
    proc::Process& process();              ///< Process object for cross-section computation/events generation
    const proc::Process& process() const;  ///< Process object for cross-section computation/events generation
    std::string processName() const;       ///< Name of the process considered
    void clearProcess();                   ///< Remove the process pointer
    void setProcess(std::unique_ptr<proc::Process>);  ///< Set a process configuration
    void setProcess(proc::Process*);                  ///< Set a process configuration

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

      inline void setTargetLuminosity(double lipb) { target_lumi_ = lipb; }  ///< Set target luminosity, in pb^-1
      inline double targetLuminosity() const { return target_lumi_; }        ///< Target luminosity to reach, in pb^-1
      inline void setMaxGen(size_t max_gen) { max_gen_ = max_gen; }  ///< Set the maximal number of events to generate
      inline size_t maxGen() const { return max_gen_; }              ///< Maximal number of events to generate
      inline bool enabled() const { return max_gen_ > 0; }           ///< Are we generating events?
      inline void setPrintEvery(size_t every) { gen_print_every_ = every; }  ///< Set the events display periodicity
      inline size_t printEvery() const { return gen_print_every_; }          ///< Periodicity of event displays
      inline void setSymmetrise(bool sym) { symmetrise_ = sym; }   ///< Symmetrisation of the z-axis for each event?
      inline bool symmetrise() const { return symmetrise_; }       ///< Symmetrise events wrt the \f$z\f$-axis ?
      inline void setNumThreads(size_t nt) { num_threads_ = nt; }  ///< Set number of threads for event generation
      inline size_t numThreads() const { return num_threads_; }    ///< Number of threads to perform event generation
      inline void setNumPoints(size_t np) { num_points_ = np; }    ///< Set number of points to probe in each integr.bin
      inline size_t numPoints() const { return num_points_; }  ///< Number of points to "shoot" in each integration bin

    private:
      int max_gen_, gen_print_every_;
      double target_lumi_;
      bool symmetrise_;
      int num_threads_, num_points_;
    };
    inline Generation& generation() { return generation_; }              ///< Event generation parameters
    inline const Generation& generation() const { return generation_; }  ///< Event generation parameters

    //----- event modification (e.g. hadronisation, decay) algorithm

    EventModifier& eventModifier(size_t) const;  ///< Event modification algorithm
    /// List of event modification algos
    inline EventModifiersSequence& eventModifiersSequence() { return evt_modifiers_; }
    /// List of event modification algos
    inline const EventModifiersSequence& eventModifiersSequence() const { return evt_modifiers_; }
    void clearEventModifiersSequence();                ///< Remove all event modifiers from sequence
    void addModifier(std::unique_ptr<EventModifier>);  ///< Add a new event modification algorithm to the sequence
    void addModifier(EventModifier*);                  ///< Add a new event modification algorithm to the sequence

    //----- event output algorithms

    EventExporter& eventExporter(size_t) const;  ///< Output module
    /// List of event output modules
    inline EventExportersSequence& eventExportersSequence() { return evt_exporters_; }
    /// List of event output modules
    inline const EventExportersSequence& eventExportersSequence() const { return evt_exporters_; }
    void clearEventExportersSequence();                     ///< Remove all output modules from sequence
    void addEventExporter(std::unique_ptr<EventExporter>);  ///< Set a new output module definition
    void addEventExporter(EventExporter*);                  ///< Set the pointer to an output module

    //----- taming functions

    /// List of all taming functions definitions
    inline const TamingFunctionsSequence& tamingFunctions() const { return taming_functions_; }
    void addTamingFunction(std::unique_ptr<utils::Functional>);  ///< Set a new taming function definition

    //----- run operations

    void prepareRun();  ///< Reset total generation time and number of events generated for this run, prepare kinematics

    /// Add a new timing into the total generation time
    /// \param[in] generation_time Time to add (in seconds)
    void addGenerationTime(double generation_time);
    inline double totalGenerationTime() const { return total_gen_time_; }  ///< Total generation time in s for this run
    inline unsigned int numGeneratedEvents() const { return num_gen_events_; }  ///< Number of events generated in run

  private:
    std::unique_ptr<proc::Process> process_;    ///< Physics process held by these parameters
    EventModifiersSequence evt_modifiers_;      ///< Collection of event modification algorithms to be applied
    EventExportersSequence evt_exporters_;      ///< Collection of event output modules to be applied
    TamingFunctionsSequence taming_functions_;  ///< Functions to be used to account for rescattering corrections
    double total_gen_time_{0.};                 ///< Total generation time (in seconds)
    unsigned long num_gen_events_{0ul};         ///< Number of events already generated
    ParametersList integrator_;                 ///< Integrator parameters
    Generation generation_;                     ///< Events generation parameters
    std::unique_ptr<utils::TimeKeeper> timer_;    ///< Collection of stopwatches for timing
  };
}  // namespace cepgen

#endif
