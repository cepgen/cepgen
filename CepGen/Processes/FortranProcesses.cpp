#include "CepGen/Processes/FortranProcesses.h"
#include "CepGen/Core/Exception.h"

//=================================================================================================
// LIST ALL FORTRAN SUBROUTINES HERE
//=================================================================================================

extern "C"
{
  extern void nucl_to_ff_( double& );
}

//=================================================================================================
BEGIN_FORTRAN_PROCESSES_ENUM // DO NOT REMOVE ME ==================================================
//=================================================================================================

REGISTER_FORTRAN_PROCESS( "patoll", nucl_to_ff_, "pA ↝ (g/ɣ)ɣ → l⁺l¯" )
REGISTER_FORTRAN_PROCESS( "patoff", nucl_to_ff_, "pA ↝ (g/ɣ)ɣ → f⁺f¯" )
REGISTER_FORTRAN_PROCESS( "aptoll", nucl_to_ff_, "Ap ↝ ɣ(g/ɣ) → l⁺l¯" )
REGISTER_FORTRAN_PROCESS( "aptoff", nucl_to_ff_, "Ap ↝ ɣ(g/ɣ) → f⁺f¯" )
REGISTER_FORTRAN_PROCESS( "aatoll", nucl_to_ff_, "AA ↝ ɣɣ → l⁺l¯" )
REGISTER_FORTRAN_PROCESS( "aatoff", nucl_to_ff_, "AA ↝ ɣɣ → f⁺f¯" )

//=================================================================================================
// DO NOT EDIT BELOW THESE LINES ==================================================================
//=================================================================================================
std::ostringstream os;
for ( auto& proc : Process::FortranProcessesHandler::get().list() )
  os << "\n *) " << proc.name << ": " << proc.description;
CG_DEBUG( "FortranProcesses" )
  << "List of Fortran processes registered:"
  << os.str();

END_FORTRAN_PROCESSES_ENUM

