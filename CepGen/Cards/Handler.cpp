#include "CepGen/Cards/PythonHandler.h"
#include "CepGen/Cards/LpairHandler.h"
#include "CepGen/Core/Exception.h"

namespace cepgen
{
  namespace card
  {
    Parameters&
    Handler::parse( const char* filename )
    {
      const std::string extension = getExtension( filename );
      if ( extension == "card" ) {
        static LpairHandler hnd( filename );
        return hnd.parameters();
      }
#ifdef PYTHON
      else if ( extension == "py" ) {
        static PythonHandler hnd( filename );
        return hnd.parameters();
      }
#endif
      throw CG_FATAL( "Cards:handler" )
        << "Failed to determine the steering card type for \"" << filename << "\"!";
    }
  }
}
