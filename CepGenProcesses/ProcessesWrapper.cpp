#include "CepGen/Modules/ProcessesFactory.h"
#include "CepGen/Processes/FortranKTProcess.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

#include "CepGen/Event/Event.h"

//=============================================================================
// START BY LISTING ALL FORTRAN FUNCTIONS
// usage:
//  DECLARE_FORTRAN_FUNCTION( function_name )
// with the Fortran function name written in lowercase
// (no trailing '_' necessary)
//=============================================================================

DECLARE_FORTRAN_FUNCTION( nucl_to_ff )

//=============================================================================
// START THE MAPPING name -> Fortran matrix element evaluation function
// usage:
//  REGISTER_FORTRAN_PROCESS( name, function_name, "description )
//=============================================================================

REGISTER_FORTRAN_PROCESS( pptoff_f77, "(p/A)(p/A) ↝ (g/ɣ)ɣ → f⁺f¯", nucl_to_ff )
