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

#ifndef CepGen_Core_EventModifier_h
#define CepGen_Core_EventModifier_h

#include <string>
#include <vector>

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Event;
  class Parameters;
  class ParametersList;
  /// Class template to interface (external/internal) events modification algorithms
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date July 2019
  class EventModifier : public NamedModule<std::string> {
  public:
    /// Default constructor for an undefined modifier
    /// \param[in] params User-controlled steering parameters for this module
    explicit EventModifier(const ParametersList& params);
    /// Virtual destructor
    virtual ~EventModifier() = default;

    static ParametersDescription parametersDescription();

    /// Set all runtime parameters steering this module
    virtual void setRuntimeParameters(const Parameters& params) { rt_params_ = &params; }
    /// \brief Specify a random numbers generator seed for the external module
    /// \param[in] seed A RNG seed
    void setSeed(long long seed) { seed_ = seed; }

    /// Parse a configuration string
    virtual void readString(const char*) {}
    /// Parse a configuration string
    virtual void readString(const std::string& param) { readString(param.c_str()); }
    /// Parse a list of configuration strings
    virtual void readStrings(const std::vector<std::string>& params);

    /// Initialise the event modifier before its running
    virtual void init() = 0;
    /** \brief Modify a full event
       * \param[inout] ev Input/output event
       * \param[inout] weight Event weight after modification
       * \param[in] full Perform the full state modification
       * \return Boolean stating whether or not the modification occurred successfully
       */
    virtual bool run(Event& ev, double& weight, bool full) = 0;
    /// Specify the process cross section and uncertainty, in pb
    virtual void setCrossSection(double, double) {}

  protected:
    /// Random numbers generator seed fed to the algorithm
    long long seed_{0ll};
    /// Maximal number of trials for the algorithm
    unsigned short max_trials_{1};
    /// List of runtime parameters steering this module
    const Parameters* rt_params_{nullptr};  // not owning
  };
}  // namespace cepgen

#endif
