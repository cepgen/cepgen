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

      const RunParameters* runParameters() const { return rt_params_; }  ///< Parsed list of runtime parameters
      RunParameters* runParameters() { return rt_params_; }              ///< Parsed list of runtime parameters

      virtual void pack(const RunParameters*) {}  ///< Specify runtime parameters to the handler

      virtual RunParameters* parseString(const std::string&, RunParameters* params) { return params; }
      /// Retrieve a configuration from a parsed steering card
      virtual RunParameters* parseFile(const std::string&, RunParameters* params) { return params; }
      static RunParameters* parseString(const std::string&);  ///< Build a configuration from a steering card
      static RunParameters* parseFile(const std::string&);    ///< Build a configuration from a steering card

      virtual void write(const std::string&) const {}  ///< Write the current configuration into a steering card
      static void write(const RunParameters*, const std::string&);  ///< Write a steering card from a configuration

    protected:
      const std::string filename_;  ///< Input filename
      RunParameters* rt_params_;    ///< List of parameters parsed from a card handler
    };
  }  // namespace card
}  // namespace cepgen

#endif
