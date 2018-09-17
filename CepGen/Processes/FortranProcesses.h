#ifndef CepGen_Processes_FortranProcesses_h
#define CepGen_Processes_FortranProcesses_h

#include "CepGen/Processes/FortranKTProcess.h"

namespace CepGen
{
  namespace Process
  {
    struct FortranProcess
    {
      const char* name;
      void ( *method )( double& );
      const char* description;
    };
    class FortranProcessesHandler
    {
      public:
        static FortranProcessesHandler& get() {
          static FortranProcessesHandler fph;
          return fph;
        }
        void add( const FortranProcess& proc ) {
          processes_.emplace_back( proc );
        }
        const std::vector<FortranProcess>& list() const { return processes_; }
        FortranProcessesHandler( const FortranProcessesHandler& ) = delete;

      private:
        explicit FortranProcessesHandler() {}
        std::vector<FortranProcess> processes_;
    };

    void generateFortranProcesses();
  }
}

#define BEGIN_FORTRAN_PROCESSES_ENUM namespace CepGen { namespace Process { void generateFortranProcesses() {
#define REGISTER_FORTRAN_PROCESS( name, method, description ) \
  CepGen::Process::FortranProcessesHandler::get().add( CepGen::Process::FortranProcess{ name, method, description } );
#define END_FORTRAN_PROCESSES_ENUM }}}

#endif

