#ifndef CepGen_Processes_ProcessesHandler_h
#define CepGen_Processes_ProcessesHandler_h

#include "CepGen/Processes/GenericProcess.h"
#include <unordered_map>
#include <memory>


#define BUILDERNM(obj) obj ## Builder
#define REGISTER_PROCESS( name, obj ) \
  class BUILDERNM( obj ) { \
    public: BUILDERNM( obj )() { CepGen::ProcessesHandler::get().registerProcess( name, new obj ); } \
  }; \
  static BUILDERNM(obj) g ## obj;

//virtual PPS::GeometryComponent* Build(std::string name) { return new obj(name); }

namespace CepGen
{
  //namespace Process { class GenericProcess; }
  class ProcessesHandler
  {
    public:
      static ProcessesHandler& get();
      ~ProcessesHandler() = default;

      void registerProcess( const std::string& name, const CepGen::Process::GenericProcess* );
      std::unique_ptr<CepGen::Process::GenericProcess> build( const std::string& name ) const;
      void dump() const;

    private:
      explicit ProcessesHandler() = default;
      std::unordered_map<std::string, std::unique_ptr<const Process::GenericProcess> > map_;

    public:
      ProcessesHandler( const ProcessesHandler& ) = delete;
      void operator=( const ProcessesHandler& ) = delete;
  };
}

#endif

