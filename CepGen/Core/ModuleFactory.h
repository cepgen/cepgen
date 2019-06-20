#ifndef CepGen_Core_ModuleFactory_h
#define CepGen_Core_ModuleFactory_h

#include "CepGen/Core/ParametersList.h"

#include <unordered_map>
#include <vector>
#include <memory>
#include <sstream>

#define BUILDERNM( obj ) obj ## Builder
#define STRINGIFY( name ) #name

namespace cepgen
{
  /// A generic factory to build modules
  /// \tparam T Base class to build
  /// \tparam I Indexing variable type
  template<typename T, typename I=std::string>
  class ModuleFactory
  {
    public:
      /// Retrieve a unique instance of this factory
      static ModuleFactory& get() {
        static ModuleFactory<T,I> instance;
        return instance;
      }
      /// Default destructor
      ~ModuleFactory() = default;
      /// Register a named module in the database
      /// \tparam U Class to register (inherited from T base class)
      template <typename U> void registerModule( const I& name, const ParametersList& def_params = ParametersList() ) {
        static_assert( std::is_base_of<T,U>::value, "Failed to register an object with improper inheritance into the factory." );
        map_[name] = &create<U>;
        params_map_[name] = def_params;
      }
      /// Build one instance of a named module
      /// \param[in] name Module name to retrieve
      /// \param[in] params List of parameters to be invoked by the constructor
      std::unique_ptr<T> build( const I& name, ParametersList params = ParametersList() ) const {
        if ( name == I() || map_.count( name ) == 0 ) {
          std::ostringstream oss;
          oss << __PRETTY_FUNCTION__ << " Failed to retrieve a module with index \"" << name << "\" from factory!";
          throw std::runtime_error( oss.str() );
        }
        if ( params_map_.count( name ) > 0 )
          params += params_map_.at( name );
        return map_.at( name )( params );
      }
      /// Build one instance of a named module
      /// \param[in] params List of parameters to be invoked by the constructor
      std::unique_ptr<T> build( ParametersList params = ParametersList() ) const {
        if ( params.has<I>( KEY ) ) {
          const I& idx = params.get<I>( KEY );
          if ( map_.count( idx ) == 0 ) {
            std::ostringstream oss;
            oss << __PRETTY_FUNCTION__ << " Failed to retrieve a module with index \"" << idx << "\" from factory!";
            throw std::runtime_error( oss.str() );
          }
          if ( params_map_.count( idx ) > 0 )
            params += params_map_.at( idx );
          return map_.at( idx )( params );
        }
        else
          throw std::runtime_error( std::string( __PRETTY_FUNCTION__ )+" Failed to retrieve an indexing key from parameters to build from factory!" );
      }
      /// List of modules registred in the database
      std::vector<I> modules() const {
        std::vector<I> out;
        for ( const auto& p : map_ )
          out.emplace_back( p.first );
        return out;
      }

      static constexpr const char* KEY = "id";

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
      std::unordered_map<I,ModCreate> map_;
      /// Database of default parameters associated to modules
      std::unordered_map<I,ParametersList> params_map_;

    public:
      ModuleFactory( const ModuleFactory& ) = delete;
      void operator=( const ModuleFactory& ) = delete;
  };
}

#endif

