#include "CepGen/Processes/ProcessesHandler.h"

//=============================================================================
// START BY LISTING ALL FORTRAN FUNCTIONS                                    //
// usage:                                                                    //
//  DECLARE_FORTRAN_FUNCTION( function_name )                                //
// with the Fortran function name written in lowercase                       //
// (no trailing '_' necessary)                                               //
//=============================================================================

DECLARE_FORTRAN_FUNCTION( nucl_to_ff )

//=============================================================================
// START THE MAPPING name -> Fortran matrix element evaluation function      //
// usage:                                                                    //
//  REGISTER_FORTRAN_PROCESS( name, function_name, "description )            //
//=============================================================================

REGISTER_FORTRAN_PROCESS( patoff, nucl_to_ff, "pA ↝ (g/ɣ)ɣ → f⁺f¯" )
REGISTER_FORTRAN_PROCESS( aptoff, nucl_to_ff, "Ap ↝ ɣ(g/ɣ) → f⁺f¯" )
REGISTER_FORTRAN_PROCESS( aatoff, nucl_to_ff, "AA ↝ ɣɣ → f⁺f¯" )

