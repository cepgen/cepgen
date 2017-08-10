#include "CepGen/Cards/Handler.h"

namespace CepGen
{
  namespace Cards
  {
    template<Type T>
    Handler<T>::Handler( const char* file )
    {
      if ( strlen( file ) != 0 ) parse( file );
    }

    template<Type T>
    Handler<T>::~Handler()
    {}
  }
}

