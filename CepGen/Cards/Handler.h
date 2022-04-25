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

#ifndef CepGen_Cards_Handler_h
#define CepGen_Cards_Handler_h

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Parameters;
  class ParametersList;
  /// Location for all steering card parsers/writers
  namespace card {
    /// Base steering card module
    class Handler : public NamedModule<std::string> {
    public:
      /// Build a configuration from an external steering card
      explicit Handler(const ParametersList&);
      virtual ~Handler() = default;

      static ParametersDescription description();

      /// Get the list of runtime parameters as parsed
      const Parameters* runtimeParameters() const { return rt_params_; }
      /// Get the list of runtime parameters as parsed
      Parameters* runtimeParameters() { return rt_params_; }
      /// Specify runtime parameters to the handler
      virtual void pack(const Parameters*){};
      /// Retrieve a configuration from a parsed steering card
      virtual Parameters* parse(const std::string&, Parameters* params) { return params; }
      /// Build a configuration from a steering card
      static Parameters* parse(const std::string&);
      /// Write the current configuration into a steering card
      virtual void write(const std::string&) const {}
      /// Write a steering card from a configuration
      static void write(const Parameters*, const std::string&);

    protected:
      /// Input filename
      const std::string filename_;
      /// List of parameters parsed from a card handler
      Parameters* rt_params_;
    };
  }  // namespace card
}  // namespace cepgen

#endif
