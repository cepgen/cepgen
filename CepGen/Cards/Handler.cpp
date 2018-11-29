#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Cards/LpairHandler.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace card
  {
    Parameters
    Handler::parse( const char* filename )
    {
      const std::string extension = getExtension( filename );
      if ( extension == "card" )
        return LpairHandler( filename ).parameters();
#ifdef PYTHON
      else if ( extension == "py" )
        return PythonHandler( filename ).parameters();
#endif
      throw CG_FATAL( "Cards:handler" )
        << "Failed to determine the steering card type for \"" << filename << "\"!";
    }
  }
}
