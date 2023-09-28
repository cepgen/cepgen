/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2016-2023  Laurent Forthomme
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

#ifndef CepGen_EventFilter_EventExporter_h
#define CepGen_EventFilter_EventExporter_h

#include <string>

#include "CepGen/EventFilter/EventHandler.h"

namespace cepgen {
  class Event;
  class Value;
  /**
   * \brief Output format handler for events export
   * \author Laurent Forthomme <laurent.forthomme@cern.ch>
   * \date Sep 2016
   */
  class EventExporter : public EventHandler {
  public:
    explicit EventExporter(const ParametersList&);

    /// Specify the process cross section and uncertainty, in pb
    virtual void setCrossSection(const Value&) {}
    /// Set the event number
    void setEventNumber(unsigned long long ev_id) { event_num_ = ev_id; }

    /// Writer operator
    virtual void operator<<(const Event&) = 0;

  protected:
    /// Print a banner containing all runtime parameters information
    std::string banner(const std::string& prep = "") const;
    /// Event index
    unsigned long long event_num_{0ull};
  };
}  // namespace cepgen

#endif
