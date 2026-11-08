/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023-2025  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/RunParameters.h"
#include "CepGen/EventFilter/EventHandler.h"

using namespace cepgen;

EventHandler::EventHandler(const ParametersList& params) : NamedModule(params) {}

EventHandler::~EventHandler() { CG_DEBUG("EventHandler") << "Destructor called for '" << name_ << "' event handler."; }

ParametersDescription EventHandler::description() {
  auto desc = ParametersDescription();
  desc.setDescription("Unnamed event handler");
  return desc;
}

void EventHandler::initialise(const RunParameters& params) {
  if (initialised_)
    CG_WARNING("EventHandler:initialise") << "Event handler '" << name_ << "' was already initialised.";
  run_params_ = &params;
  initialise();
  initialised_ = true;
}

const RunParameters& EventHandler::runParameters() const {
  if (!run_params_)
    throw CG_FATAL("EventHandler:runParameters") << "Run parameters not yet initialised.";
  return *run_params_;
}

void* EventHandler::enginePtr() {
  throw CG_FATAL("EventHandler:enginePtr")
      << "No engine object declared for event handler with name '" << name_ << "'.";
}
