#include "CepGen/Event/Event.h"
#include "CepGen/Modules/ProcessFactory.h"
#include "CepGen/Processes/FortranKTProcess.h"
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

REGISTER_FORTRAN_PROCESS(pptoff_f77, "(p/A)(p/A) ↝ (g/ɣ)ɣ → f⁺f¯", nucl_to_ff)
