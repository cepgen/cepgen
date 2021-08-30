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

#ifndef CepGen_Modules_ExportModule_h
#define CepGen_Modules_ExportModule_h

#include <iosfwd>
#include <string>

#include "CepGen/Modules/NamedModule.h"

namespace cepgen {
  class Event;
  class Parameters;
  /// Location for all output generators
  namespace io {
    /**
     * \brief Output format handler for events export
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Sep 2016
     */
    class ExportModule : public NamedModule<std::string> {
    public:
      /// Class constructor
      /// \param[in] params User-controlled steering parameters for this module
      explicit ExportModule(const ParametersList& params);
      virtual ~ExportModule();

      /// Set the process cross section and its associated error
      virtual void setCrossSection(double /*cross_section*/, double /*err_cross_section*/) {}
      /// Set the event number
      void setEventNumber(const unsigned int& ev_id) { event_num_ = ev_id; }

      /// Initialise the handler and its inner parameterisation
      virtual void initialise(const Parameters&) = 0;
      /// Writer operator
      virtual void operator<<(const Event&) = 0;

    protected:
      /// Print a banner containing all runtime parameters information
      static std::string banner(const Parameters&, const std::string& prep = "");
      /// Event index
      unsigned long long event_num_{0ull};
    };
  }  // namespace io
}  // namespace cepgen

#endif
