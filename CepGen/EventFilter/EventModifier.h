/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2019-2025  Laurent Forthomme
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

#ifndef CepGen_EventFilter_EventModifier_h
#define CepGen_EventFilter_EventModifier_h

#include "CepGen/EventFilter/EventHandler.h"

namespace cepgen {
  class Value;
  /// Class template to interface (external/internal) events modification algorithms
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date July 2019
  class EventModifier : public EventHandler {
  public:
    explicit EventModifier(const ParametersList&);  ///< Default constructor for an undefined modifier

    static ParametersDescription description();

    /// Specify a random numbers generator seed for the external module
    /// \param[in] seed An RNG seed
    void setSeed(long long seed) { seed_ = seed; }

    inline virtual void readString(const std::string&) {}       ///< Parse a configuration string
    virtual void readStrings(const std::vector<std::string>&);  ///< Parse a list of configuration strings

    /// Modify an event
    /// \param[inout] ev Input/output event
    /// \param[inout] weight Event weight after modification
    /// \param[in] fast Run a faster version of the algorithm (whenever available)
    /// \return Boolean stating whether the modification occurred successfully
    virtual bool run(Event& ev, double& weight, bool fast = false) = 0;
    inline virtual void setCrossSection(const Value&) {}  ///< Specify the cross-section value, in pb

  protected:
    long long seed_{0ll};           ///< Random numbers generator seed fed to the algorithm
    unsigned short max_trials_{1};  ///< Maximal trials for the algorithm
  };
}  // namespace cepgen

#endif
