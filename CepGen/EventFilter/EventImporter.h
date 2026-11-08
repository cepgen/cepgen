/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2022-2025  Laurent Forthomme
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

#include "CepGen/EventFilter/EventHandler.h"
#include "CepGen/Utils/Value.h"

namespace cepgen {
  /// Base event importer module
  /// \author Laurent Forthomme <laurent.forthomme@cern.ch>
  /// \date Dec 2022
  class EventImporter : public EventHandler {
  public:
    using EventHandler::EventHandler;

    static ParametersDescription description() {
      auto desc = EventHandler::description();
      desc.setDescription("Unnamed event importer");
      return desc;
    }

    virtual bool operator>>(Event&) = 0;                          ///< Read the next event
    const Value& crossSection() const { return cross_section_; }  ///< Process cross-section and uncertainty, in pb

  protected:
    /// Specify the process cross-section and uncertainty, in pb
    void setCrossSection(const Value& cross_section) { cross_section_ = cross_section; }

  private:
    Value cross_section_;
  };
}  // namespace cepgen

#endif
