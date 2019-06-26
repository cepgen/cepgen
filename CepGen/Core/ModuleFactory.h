#ifndef CepGen_Core_ModuleFactory_h
#define CepGen_Core_ModuleFactory_h

#include <unordered_map>
#include <vector>
#include <memory>
#include <iosfwd>

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name

namespace cepgen
{
  class ParametersList;
  /// A generic factory to build modules
  /// \tparam T Base class to build
  template<typename T>
  class ModuleFactory
  {
    public:
      /// Retrieve a unique instance of this factory
      static ModuleFactory& get() {
        static ModuleFactory<T> instance;
        return instance;
      }
      /// Default destructor
      ~ModuleFactory() = default;
      /// Register a named module in the database
      /// \tparam U Class to register (inherited from T base class)
      template <typename U> void registerModule( const std::string& name ) {
        static_assert( std::is_base_of<T,U>::value, "Failed to register a class of non-base class template." );
        map_[name] = &create<U>;
      }
      /// Build one instance of a named module
      /// \param[in] name Module name to retrieve
      /// \param[in] params List of parameters to be invoked by the constructor
      std::unique_ptr<T> build( const std::string& name, const ParametersList& params = ParametersList() ) const {
        if ( name.empty() || map_.count( name ) == 0 )
          throw std::runtime_error( std::string( __PRETTY_FUNCTION__ )+" Failed to retrieve a module with name \""+name+"\" from factory!" );
        return map_.at( name )( params );
      }
      /// List of modules registred in the database
      std::vector<std::string> modules() const {
        std::vector<std::string> out;
        for ( const auto& p : map_ )
          out.emplace_back( p.first );
        return out;
      }

    private:
      explicit ModuleFactory() = default;
      /// Construct a module with its parameters set
      template<typename U> static std::unique_ptr<T> create( const ParametersList& params ) {
        return std::unique_ptr<T>( new U( params ) );
      }

    protected:
      /// Constructor type for a module
      typedef std::unique_ptr<T> (*ModCreate)( const ParametersList& );
      /// Database of modules handled by this instance
      std::unordered_map<std::string,ModCreate> map_;

    public:
      ModuleFactory( const ModuleFactory& ) = delete;
      void operator=( const ModuleFactory& ) = delete;
  };
}

#endif

