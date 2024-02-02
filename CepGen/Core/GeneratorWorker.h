/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#include <memory>
#include <vector>

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Event/Event.h"

namespace cepgen {
  class Integrator;
  class RunParameters;
  class ProcessIntegrand;
  namespace proc {
    class Process;
  }
  /// Monte-Carlo generator instance
  class GeneratorWorker : public SteeredObject<GeneratorWorker> {
  public:
    /// Book the memory slots and structures for the generator
    explicit GeneratorWorker(const ParametersList&);
    virtual ~GeneratorWorker();

    static ParametersDescription description();

    /// Specify the runtime parameters
    void setRunParameters(const RunParameters*);
    /// Specify the integrator instance handled by the mother generator
    void setIntegrator(const Integrator* integ);
    /// Launch the event generation
    /// \param[in] num_events Number of events to generate
    /// \param[in] callback The callback function applied on every event generated
    void generate(size_t num_events, const std::function<void(const proc::Process&)>&);
    /// Function evaluator
    ProcessIntegrand& integrand() { return *integrand_; }

    /// Initialise the generation parameters
    virtual void initialise() = 0;
    /// Generate a single event
    virtual bool next() = 0;

  protected:
    /// Store the event in the output file
    /// \param[in] callback The callback function for every event generated
    /// \return A boolean stating whether or not the event was successfully saved
    bool storeEvent();

    /// Pointer to the mother-handled integrator instance
    /// \note NOT owning
    const Integrator* integrator_{nullptr};
    /// Steering parameters for the event generation
    /// \note NOT owning
    const RunParameters* params_{nullptr};
    /// Local event weight evaluator
    std::unique_ptr<ProcessIntegrand> integrand_;
    /// Callback function on process for each new event
    std::function<void(const proc::Process&)> callback_proc_{nullptr};
  };
}  // namespace cepgen

#endif
