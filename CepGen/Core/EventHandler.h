/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#ifndef CepGen_Core_EventHandler_h
#define CepGen_Core_EventHandler_h

#include <iosfwd>
#include <string>
#include <vector>

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Event;
  class Parameters;
  /// Class template for modules interacting with events
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Jan 2023
  class EventHandler : public NamedModule<std::string> {
  public:
    explicit EventHandler(const ParametersList&);
    virtual ~EventHandler();

    static ParametersDescription description();

    /// Initialise the handler and its inner parameterisation
    void initialise(const Parameters&);
    /// List of run parameters
    const Parameters& runParameters() const;

    /// Retrieve the engine object
    template <typename T>
    T* engine() {
      return static_cast<T*>(enginePtr());
    }

  protected:
    virtual void initialise() = 0;
    /// Engine object
    virtual void* enginePtr();

  private:
    const Parameters* run_params_{nullptr};  // NOT owning
    bool initialised_{false};
  };
}  // namespace cepgen

#endif
