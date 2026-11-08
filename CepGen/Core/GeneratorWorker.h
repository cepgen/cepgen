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

#include <memory>

#include "CepGen/Event/Event.h"

namespace cepgen {
  class Integrator;
  class RunParameters;
  class ProcessIntegrand;
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

    void setRunParameters(const RunParameters*);  ///< Specify the runtime parameters
    void setIntegrator(const Integrator*);        ///< Specify the integrator instance handled by the mother generator

    /// Launch the event generation
    /// \param[in] num_events Events multiplicity to generate
    /// \param[in] callback The callback function applied on every event generated
    void generate(size_t num_events, const std::function<void(const proc::Process&)>& callback);

    inline ProcessIntegrand& integrand() const { return *integrand_; }  ///< Function evaluator

    virtual void initialise() = 0;  ///< Initialise the generation parameters
    virtual bool next() = 0;        ///< Generate a single event

  protected:
    /// Store the event in the output file
    /// \return A boolean stating whether the event was successfully saved
    bool storeEvent() const;

    // NOT owned
    const Integrator* integrator_{nullptr};     ///< Pointer to the mother-handled integrator instance
    const RunParameters* run_params_{nullptr};  ///< Steering parameters for the event generation

    std::unique_ptr<ProcessIntegrand> integrand_;                       ///< Local event weight evaluator
    std::function<void(const proc::Process&)> callback_proc_{nullptr};  ///< Callback function for each new event
  };
}  // namespace cepgen

#endif
