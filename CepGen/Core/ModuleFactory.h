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
      static ModuleFactory& get();
      ~ModuleFactory() = default;

      void registerModule( const std::string& name, const T* );
      std::unique_ptr<T> build( const std::string& name, const ParametersList& ) const;
      void dump() const;

    private:
      explicit ModuleFactory() = default;
      std::unordered_map<std::string, std::unique_ptr<const T> > map_;

    public:
      ModuleFactory( const ModuleFactory& ) = delete;
      void operator=( const ModuleFactory& ) = delete;
  };
}

#endif
