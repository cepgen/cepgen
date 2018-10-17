#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Processes/ProcessesHandler.h"

#include <sstream>

namespace cepgen
{
  template<typename T>
  ModuleFactory<T>&
  ModuleFactory<T>::get()
  {
    static ModuleFactory<T> instance;
    return instance;
  }

  template<typename T> void
  ModuleFactory<T>::registerModule( const std::string& name, const T* proc )
  {
    map_[name].reset( proc );
    CG_DEBUG( "ModuleFactory" ) << "Module name \"" << name << "\" registered in database.";
  }

  template<typename T> void
  ModuleFactory<T>::dump() const
  {
    std::ostringstream oss;
    for ( const auto& p : map_ )
      oss << " '" << p.first << "'";
    CG_INFO( "ModuleFactory:dump" )
      << "List of process(es) handled in the database:" << oss.str();
  }

  template<typename T> std::unique_ptr<T>
  ModuleFactory<T>::build( const std::string& name, const ParametersList& params ) const
  {
    if ( map_.count( name ) == 0 )
      throw CG_FATAL( "ModuleFactory:build" )
        << "Failed to retrieve a process with name \"" << name << "\"!";
    return map_.at( name )->clone( params );
  }
}
