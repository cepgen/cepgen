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

#ifndef CepGenAddOns_PythonWrapper_Environment_h
#define CepGenAddOns_PythonWrapper_Environment_h

#include "CepGen/Core/SteeredObject.h"
#include "CepGenAddOns/PythonWrapper/PythonTypes.h"

namespace cepgen {
  namespace python {
    class Environment : SteeredObject<Environment> {
    public:
      /// Initialise the python environment
      explicit Environment(const ParametersList&);
      /// Finalise the python environment
      ~Environment();

      static ParametersDescription description();

      /// Set the name of the Python program
      void setProgramName(const std::string&);
      /// Is the python environment already initialised?
      bool initialised();

    private:
#if PY_VERSION_HEX >= 0x03080000
      PyConfig config_;
#endif
    };
  }  // namespace python
}  // namespace cepgen

#endif
