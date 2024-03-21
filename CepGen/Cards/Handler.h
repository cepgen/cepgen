/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
  class RunParameters;
  /// Location for all steering card parsers/writers
  namespace card {
    /// Base steering card module
    class Handler : public NamedModule<Handler, std::string> {
    public:
      explicit Handler(const ParametersList&);  ///< Build a configuration from an external steering card
      virtual ~Handler() = default;

      static ParametersDescription description();

      /// Retrieve a configuration from a parsed steering string
      inline virtual Handler& parseCommands(const std::vector<std::string>&) { return *this; }
      /// Retrieve a configuration from a parsed steering card
      inline virtual Handler& parseFile(const std::string&) { return *this; }

      inline const RunParameters* runParameters() const { return rt_params_; }  ///< Parsed list of runtime parameters
      inline RunParameters* runParameters() { return rt_params_; }              ///< Parsed list of runtime parameters
      /// Specify runtime parameters to the handler
      inline virtual Handler& setRunParameters(const RunParameters* params) {
        rt_params_ = const_cast<RunParameters*>(params);
        return *this;
      }
      /// Write a steering card from a configuration
      inline virtual void write(const std::string&) const {}

    protected:
      const std::string filename_;  ///< Input filename

    private:
      RunParameters* rt_params_{nullptr};  ///< List of parameters parsed from a card handler
    };
  }  // namespace card
}  // namespace cepgen

#endif
