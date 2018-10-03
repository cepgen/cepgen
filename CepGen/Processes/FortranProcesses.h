#ifndef CepGen_Processes_FortranProcesses_h
#define CepGen_Processes_FortranProcesses_h

#include "CepGen/Processes/FortranKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    /// A Fortran process handler
    struct FortranProcess
    {
      const char* name; ///< CepGen-readable process name
      void ( *method )( double& ); ///< Pointer to the weight computation functional
      const char* description; ///< Human-readable process description
    };
    /// Fortran processes collector
    class FortranProcessesHandler
    {
      public:
        /// Static collector retrieval method
        static FortranProcessesHandler& get() {
          static FortranProcessesHandler fph;
          return fph;
        }
        /// Register a Fortran process into the collector
        void add( const FortranProcess& proc ) {
          processes_.emplace_back( proc );
        }
        /// Get a list of processes handled by this collector
        const std::vector<FortranProcess>& list() const { return processes_; }
        /// Generic copy-constructor
        FortranProcessesHandler( const FortranProcessesHandler& ) = delete;

      private:
        explicit FortranProcessesHandler() {}
        std::vector<FortranProcess> processes_;
    };

    void generateFortranProcesses();
  }
}
#define DECLARE_FORTRAN_SUBROUTINE( method ) \
  extern "C" { extern void method ## _( double& ); }
#define BEGIN_FORTRAN_PROCESSES_ENUM \
  namespace CepGen { namespace Process { void generateFortranProcesses() {
#define REGISTER_FORTRAN_PROCESS( name, method, description ) \
  CepGen::Process::FortranProcessesHandler::get().add( CepGen::Process::FortranProcess{ name, method ## _, description } );
#define END_FORTRAN_PROCESSES_ENUM }}}

#endif

