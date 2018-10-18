#include "CepGen/Processes/ProcessesHandler.h"

//=================================================================================================
// START BY LISTING ALL FORTRAN SUBROUTINES                                                      //
// usage:                                                                                        //
//  DECLARE_FORTRAN_SUBROUTINE( subroutine_name )                                                //
// with the Fortran subroutine name written in lowercase (no trailing '_' necessary)             //
//=================================================================================================

DECLARE_FORTRAN_SUBROUTINE( nucl_to_ff )

//=================================================================================================
// START THE MAPPING name -> Fortran SUBROUTINE                                                  //
// usage:                                                                                        //
//  REGISTER_FORTRAN_PROCESS( name, subroutine_name, "description )                              //
//=================================================================================================

REGISTER_FORTRAN_PROCESS( patoff, nucl_to_ff, "pA ↝ (g/ɣ)ɣ → f⁺f¯" )
REGISTER_FORTRAN_PROCESS( aptoff, nucl_to_ff, "Ap ↝ ɣ(g/ɣ) → f⁺f¯" )
REGISTER_FORTRAN_PROCESS( aatoff, nucl_to_ff, "AA ↝ ɣɣ → f⁺f¯" )

