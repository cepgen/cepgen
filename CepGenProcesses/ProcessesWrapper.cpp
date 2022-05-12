/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Process/FortranKTProcess.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

//=============================================================================
// START BY LISTING ALL FORTRAN FUNCTIONS
// usage:
//  DECLARE_FORTRAN_FUNCTION(function_name)
// with the Fortran function name written in lowercase (no trailing '_')
//=============================================================================

DECLARE_FORTRAN_FUNCTION(nucl_to_ff)

//=============================================================================
// START THE MAPPING name -> Fortran matrix element evaluation function
// usage:
//  REGISTER_FORTRAN_PROCESS(name, "description", function_name)
//=============================================================================

REGISTER_FORTRAN_PROCESS(pptoff_f77, "(p/A)(p/A) ↝ (g/γ)γ → f⁺f¯", nucl_to_ff)
