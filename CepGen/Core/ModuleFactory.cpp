#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Core/ParametersList.h"

#include "CepGen/Core/ModuleFactory.h"

#include <sstream>

namespace cepgen
{
  template<typename T> void
  ModuleFactory<T>::dump() const
  {
    std::ostringstream oss;
    for ( const auto& p : map_ )
      oss << " '" << p.first << "'";
    CG_INFO( "ModuleFactory:dump" )
      << "List of process(es) handled in the database:" << oss.str();
  }
}
