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

#include "CepGen/Core/EventHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Parameters.h"

namespace cepgen {
  EventHandler::EventHandler(const ParametersList& params) : NamedModule(params) {}

  EventHandler::~EventHandler() {
    CG_DEBUG("EventHandler") << "Destructor called for '" << name_ << "' event handler.";
  }

  ParametersDescription EventHandler::description() {
    auto desc = ParametersDescription();
    return desc;
  }

  void EventHandler::initialise(const Parameters& params) {
    if (initialised_)
      CG_WARNING("EventHandler:initialise") << "Event handler '" << name_ << "' was already initialised.";
    rt_params_ = &params;
    initialise();
    initialised_ = true;
  }

  void* EventHandler::engine() {
    throw CG_FATAL("EventHandler") << "No engine declared for event handler with name '" << name_ << "'.";
  }
}  // namespace cepgen
