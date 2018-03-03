#ifndef CepGen_Processes_ProcessesHandler_h
#define CepGen_Processes_ProcessesHandler_h

#include <map>

namespace CepGen
{
  class GenericProcess;
  class ProcessesHandler
  {
    public:
      static ProcessesHandler& get();
      ~ProcessesHandler();

      CepGen::Process::GenericProcess* build( const char* name ) const;

    private:
      explicit ProcessesHandler();
      std::map<const char*, Process::GenericProcess*> map_;

    public:
      ProcessesHandler( const ProcessesHandler& ) = delete;
      void operator=( const ProcessesHandler& ) = delete;
  };
}

#endif

