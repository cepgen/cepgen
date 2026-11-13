/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2026  Laurent Forthomme
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

#ifndef CepGen_Process_IntegrationVariable_h
#define CepGen_Process_IntegrationVariable_h

#include <functional>
#include <string>
#include <vector>

#include "CepGen/Process/Process.h"
#include "CepGen/Utils/Limits.h"

namespace cepgen {
  /// Description of one integration variable declared by a phase space building block
  /// (single-beam parton sampler, central system phase space, ...), to be registered by the consumer process
  struct IntegrationVariable {
    std::reference_wrapper<double> value;  ///< Reference to the variable to be mapped
    proc::Process::Mapping mapping;        ///< Mapping type
    Limits range;                          ///< Integration range (in mapped coordinates)
    std::string name;                      ///< Computer-readable variable name
    std::string description;               ///< Human-readable variable description
  };

  /// Register a single integration variable in a process
  inline void defineVariable(proc::Process& process, const IntegrationVariable& var) {
    process.defineVariable(var.value.get(), var.mapping, var.range, var.name, var.description);
  }
  /// Register a list of integration variables in a process (in the list order)
  inline void defineVariables(proc::Process& process, const std::vector<IntegrationVariable>& vars) {
    for (const auto& var : vars)
      defineVariable(process, var);
  }
}  // namespace cepgen

#endif
