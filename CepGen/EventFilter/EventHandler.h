/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2024  Laurent Forthomme
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

#ifndef CepGen_EventFilter_EventHandler_h
#define CepGen_EventFilter_EventHandler_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Event;
  class RunParameters;
  /// Class template for modules interacting with events
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jan 2023
  class EventHandler : public NamedModule<EventHandler> {
  public:
    explicit EventHandler(const ParametersList&);
    ~EventHandler() override;

    static ParametersDescription description();

    void initialise(const RunParameters&);       ///< Initialise the handler and its inner parameterisation
    const RunParameters& runParameters() const;  ///< List of run parameters

    /// Retrieve the engine object
    template <typename T>
    T* engine() {
      return static_cast<T*>(enginePtr());
    }

  protected:
    virtual void initialise() = 0;
    virtual void* enginePtr();  ///< Engine object

  private:
    const RunParameters* run_params_{nullptr};  // NOT owning
    bool initialised_{false};
  };
}  // namespace cepgen

#endif
