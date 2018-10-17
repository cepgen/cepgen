#ifndef CepGen_Processes_ProcessesHandler_h
#define CepGen_Processes_ProcessesHandler_h

#include "CepGen/Processes/GenericProcess.h"
#include <unordered_map>
#include <memory>

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name
#define REGISTER_PROCESS( name, obj ) \
  struct BUILDERNM( name ) { \
    BUILDERNM( name )() { cepgen::ModuleFactory<cepgen::proc::GenericProcess>::get().registerModule( STRINGIFY( name ), new obj ); } }; \
  static BUILDERNM( name ) g ## name;

namespace cepgen
{
  class ParametersList;
  template<typename T>
  class ModuleFactory
  {
    public:
      static ModuleFactory& get();
      ~ModuleFactory() = default;

      void registerModule( const std::string& name, const T* );
      std::unique_ptr<T> build( const std::string& name, const ParametersList& ) const;
      void dump() const;

    private:
      explicit ModuleFactory() = default;
      std::unordered_map<std::string, std::unique_ptr<const T> > map_;

    public:
      ModuleFactory( const ModuleFactory& ) = delete;
      void operator=( const ModuleFactory& ) = delete;
  };
  typedef ModuleFactory<cepgen::proc::GenericProcess> ProcessesHandler;
}

#endif
