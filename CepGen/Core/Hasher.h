#ifndef CepGen_Core_Hasher_h
#define CepGen_Core_Hasher_h

#include <cstddef>
#include <functional>

namespace cepgen
{
  /// A hasher table for a given structure
  template<class T,bool>
  struct Hasher
  {
    inline size_t operator()( const T& t ) const {
      return std::hash<T>()( t );
    }
  };
  /// A hasher table for a given structure
  template<class T>
  struct Hasher<T, true>
  {
    inline size_t operator() ( const T& t ) {
      typedef typename std::underlying_type<T>::type enumType;
      return std::hash<enumType>()( static_cast<enumType>( t ) );
    }
  };
  /// A hasher table for an enumeration
  template<class T>
  struct EnumHash
  {
    inline size_t operator()( const T& t ) const {
      return Hasher<T,std::is_enum<T>::value>()( t );
    }
  };
}

#endif

