#ifndef CepGen_Core_ModuleFactory_h
#define CepGen_Core_ModuleFactory_h

#include <unordered_map>
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

      void registerModule( const std::string& name, const T* proc ) {
        map_[name].reset( proc );
      }
      virtual std::unique_ptr<T> build( const std::string& name, const ParametersList& params ) const {
        if ( map_.count( name ) == 0 )
          throw std::logic_error( "Failed to retrieve a process with name \""+name+"\"!" );
        return map_.at( name )->clone( params );
      }
      void dump() const;

    private:
      explicit ModuleFactory() = default;

    protected:
      std::unordered_map<std::string, std::unique_ptr<const T> > map_;

    public:
      ModuleFactory( const ModuleFactory& ) = delete;
      void operator=( const ModuleFactory& ) = delete;
  };
}

#endif

