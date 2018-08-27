#ifndef CepGen_Core_Hasher_h
#define CepGen_Core_Hasher_h

namespace CepGen
{
  template<class T,bool>
  struct hasher
  {
    inline size_t operator()( const T& t ) const {
      return std::hash<T>()( t );
    }
  };
  template<class T>
  struct hasher<T, true>
  {
    inline size_t operator() ( const T& t ) {
      typedef typename std::underlying_type<T>::type enumType;
      return std::hash<enumType>()( static_cast<enumType>( t ) );
    }
  };
  template<class T>
  struct EnumHash
  {
    inline size_t operator()( const T& t ) const {
      return hasher<T,std::is_enum<T>::value>()( t );
    }
  };
}

#endif

