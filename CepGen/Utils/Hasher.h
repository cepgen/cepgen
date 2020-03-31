#ifndef CepGen_Core_Hasher_h
#define CepGen_Core_Hasher_h

#include <cstddef>
#include <functional>

namespace cepgen
{
  namespace utils
  {
    /// A hasher table for a given structure
    template<class T,bool>
    struct Hasher
    {
      /// Hash a generic table
      inline size_t operator()( const T& t ) const {
        return std::hash<T>()( t );
      }
    };
    /// A hasher table for a given structure
    template<class T>
    struct Hasher<T, true>
    {
      /// Hash a structure-indexed table
      inline size_t operator() ( const T& t ) {
        typedef typename std::underlying_type<T>::type enumType;
        return std::hash<enumType>()( static_cast<enumType>( t ) );
      }
    };
    /// A hasher table for an enumeration
    template<class T>
    struct EnumHash
    {
      /// Hash an enumerator-indexed table
      inline size_t operator()( const T& t ) const {
        return Hasher<T,std::is_enum<T>::value>()( t );
      }
    };
  }
}

#endif

