#ifndef CepGen_Core_ModuleFactory_h
#define CepGen_Core_ModuleFactory_h

#include <unordered_map>
#include <vector>
#include <memory>

namespace cepgen
{
  class ParametersList;

  template<typename T>
  class ModuleFactory
  {
    public:
      static ModuleFactory& get() {
        static ModuleFactory<T> instance;
        return instance;
      }
      ~ModuleFactory() = default;

      template <typename U> void registerModule( const std::string& name ) {
        static_assert( std::is_base_of<T,U>::value, "Failed to register a class of non-base class template." );
        map_[name] = &create<U>;
      }
      std::unique_ptr<T> build( const std::string& name, const ParametersList& params ) const {
        if ( map_.count( name ) == 0 )
          throw std::runtime_error( "Failed to retrieve a module with name \""+name+"\" from factory!" );
        return map_.at( name )( params );
      }
      std::vector<std::string> modules() const {
        std::vector<std::string> out;
        for ( const auto& p : map_ )
          out.emplace_back( p.first );
        return out;
      }

    private:
      explicit ModuleFactory() = default;
      template<typename U> static std::unique_ptr<T> create( const ParametersList& params ) {
        return std::unique_ptr<T>( new U( params ) );
      }

    protected:
      typedef std::unique_ptr<T> (*ModCreate)( const ParametersList& );
      std::unordered_map<std::string,ModCreate> map_;

    public:
      ModuleFactory( const ModuleFactory& ) = delete;
      void operator=( const ModuleFactory& ) = delete;
  };
}

#endif

