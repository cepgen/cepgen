/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

#ifndef CepGen_Core_GeneratorWorker_h
#define CepGen_Core_GeneratorWorker_h

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>

#include "CepGen/Event/Event.h"

namespace cepgen {
  class Integrator;
  class RunParameters;
  class ProcessIntegrand;
  class Value;
}  // namespace cepgen
namespace cepgen::proc {
  class Process;
}  // namespace cepgen::proc

namespace cepgen {
  /// Event generator worker instance
  class GeneratorWorker : public SteeredObject<GeneratorWorker> {
  public:
    explicit GeneratorWorker(const ParametersList&);  ///< Book memory slots and structures for the generator
    ~GeneratorWorker() override;

    static ParametersDescription description();

    void setRunParameters(const RunParameters*, size_t stream = 0);  ///< Specify the stream-local runtime parameters
    void setIntegrator(const Integrator*);  ///< Specify the integrator instance handled by the mother generator

    /// Launch the event generation
    /// \param[in] num_events Events multiplicity to generate
    /// \param[in] callback The callback function applied on every event generated
    void generate(size_t num_events, const std::function<void(const proc::Process&)>& callback);

    /// Generate a fixed number of events, to be called from one of the generation threads
    /// \param[in] stop Flag requesting this worker to stop after the event being generated
    void produce(size_t num_events,
                 const std::function<void(const proc::Process&)>& callback,
                 const std::atomic<bool>& stop);

    /// Build an independent worker for a secondary generation thread, sharing this worker's preparation
    /// \param[in] stream Index of the thread (strictly positive)
    /// \param[in] cross_section Cross-section computed for this run
    /// \return A worker ready to generate events, or a null pointer if this algorithm cannot be multithreaded
    inline virtual std::unique_ptr<GeneratorWorker> fork(size_t /*stream*/, const Value& /*cross_section*/) const {
      return nullptr;
    }
    void setStoreMutex(std::mutex* store_mutex) { store_mutex_ = store_mutex; }

    inline ProcessIntegrand& integrand() const { return *integrand_; }  ///< Function evaluator

    virtual void initialise() = 0;      ///< Initialise the generation parameters
    virtual bool next() = 0;            ///< Generate a single event
    virtual bool weighted() const = 0;  ///< Are events generated cross section-weighted?

  protected:
    /// Store the event in the output file
    /// \param[in] weight Event weight, normalised to a unit mean (1 for unweighted events)
    /// \return A boolean stating whether the event was successfully saved
    bool storeEvent(double weight = 1.) const;

    // NOT owned
    const Integrator* integrator_{nullptr};     ///< Pointer to the mother-handled integrator instance
    const RunParameters* run_params_{nullptr};  ///< Steering parameters for the event generation

    std::unique_ptr<ProcessIntegrand> integrand_;                       ///< Local event weight evaluator
    std::function<void(const proc::Process&)> callback_proc_{nullptr};  ///< Callback function for each new event
    std::mutex* store_mutex_{nullptr};  ///< NOT owned: lock shared by all workers for the events storage
  };
}  // namespace cepgen

#endif
