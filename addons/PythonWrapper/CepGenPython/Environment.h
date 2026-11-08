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

#ifndef CepGenPython_Environment_h
#define CepGenPython_Environment_h

#include <Python.h>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen::python {
  class Environment : SteeredObject<Environment> {
  public:
    explicit Environment(const ParametersList&);  ///< Initialise the python environment
    ~Environment() override;                      ///< Finalise the python environment

    static ParametersDescription description();

    void setProgramName(const std::string&);  ///< Set the name of the Python program
    static bool initialised();                ///< Is the python environment already initialised?
#if PY_VERSION_HEX >= 0x03080000
    inline const PyConfig& configuration() const { return config_; }  ///< Set of configuration variables
#endif

  private:
#if PY_VERSION_HEX >= 0x03080000
    PyConfig config_;
#endif
  };
}  // namespace cepgen::python

#endif
