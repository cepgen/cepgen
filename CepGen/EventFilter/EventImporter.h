/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2023  Laurent Forthomme
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

#ifndef CepGen_EventFilter_EventImporter_h
#define CepGen_EventFilter_EventImporter_h

#include "CepGen/Event/Event.h"
#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Event;
  /// Base event importer module
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Dec 2022
  class EventImporter : public NamedModule<std::string> {
  public:
    explicit EventImporter(const ParametersList& params) : NamedModule(params) {}

    static ParametersDescription description() {
      auto desc = ParametersDescription();
      desc.setDescription("Unnamed event importer");
      return desc;
    }

    /// Output format-custom extraction operator
    template <typename T>
    Event convert(const T& in) const {
      Event out;
      convert((void*)&in, out);
      return out;
    }

    /// Read the next event
    virtual bool next(Event&) const { return false; }

  private:
    /// Output format-custom conversion algorithm
    virtual void convert(const void*, Event&) const {}
  };
}  // namespace cepgen

#endif
