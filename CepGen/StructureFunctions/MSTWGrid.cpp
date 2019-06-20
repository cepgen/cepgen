#include "CepGen/StructureFunctions/MSTWGrid.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include <fstream>

namespace mstw
{
  Grid::Grid( const cepgen::ParametersList& params ) :
    cepgen::strfun::Parameterisation( params ),
    cepgen::GridHandler<2,2>( cepgen::GridType::logarithmic )
  {
    { // file readout part
      const std::string grid_path = params_.get<std::string>( "gridPath", DEFAULT_MSTW_GRID_PATH );
      std::ifstream file( grid_path, std::ios::binary | std::ios::in );
      if ( !file.is_open() )
        throw CG_FATAL( "MSTW" ) << "Failed to load grid file \"" << grid_path << "\"!";

      file.read( reinterpret_cast<char*>( &header_ ), sizeof( header_t ) );

      // first checks on the file header

      if ( header_.magic != GOOD_MAGIC )
        throw CG_FATAL( "MSTW" ) << "Wrong magic number retrieved: " << header_.magic << ", expecting " << GOOD_MAGIC << ".";

      if ( header_.nucleon != header_t::proton )
        throw CG_FATAL( "MSTW" ) << "Only proton structure function grids can be retrieved for this purpose!";

      // retrieve all points and evaluate grid boundaries

      sfval_t val;
      while ( file.read( reinterpret_cast<char*>( &val ), sizeof( sfval_t ) ) )
        insert( { val.xbj, val.q2 }, { val.f2, val.fl } );
      file.close();
    }

    init();

    const auto& bounds = boundaries();
    CG_DEBUG( "MSTW" )
      << "MSTW@" << header_.order << " grid evaluator built "
      << "for " << header_.nucleon << " structure functions (" << header_.cl << ")\n\t"
      << "xBj in range [" << pow( 10., bounds[0].first ) << ":" << pow( 10., bounds[0].second ) << "]\n\t"
      << " Q² in range [" << pow( 10., bounds[1].first ) << ":" << pow( 10., bounds[1].second ) << "].";
  }

  std::string
  Grid::description() const
  {
    std::ostringstream os;
    const auto& bounds = boundaries();
    os << "MSTW(grid){"
       << pow( 10., bounds[0].first ) << "<xbj<" << pow( 10., bounds[0].second ) << ","
       << pow( 10., bounds[1].first ) << "<Q²/GeV²<" << pow( 10., bounds[1].second ) << "}";
    return os.str();
  }

  Grid&
  Grid::operator()( double xbj, double q2 )
  {
    const std::array<double,2> val = cepgen::GridHandler<2,2>::eval( { xbj, q2 } );
    F2 = val[0];
    FL = val[1];
    return *this;
  }

  std::ostream&
  operator<<( std::ostream& os, const Grid::sfval_t& val )
  {
    return os << cepgen::Form( "xbj = %.4f\tQ² = %.5e GeV²\tF₂ = % .6e\tFₗ = % .6e", val.xbj, val.q2, val.f2, val.fl );
  }

  std::ostream&
  operator<<( std::ostream& os, const Grid::header_t::order_t& order )
  {
    switch ( order ) {
      case Grid::header_t::lo: return os << "LO";
      case Grid::header_t::nlo: return os << "nLO";
      case Grid::header_t::nnlo: return os << "nnLO";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Grid::header_t::cl_t& cl )
  {
    switch ( cl ) {
      case Grid::header_t::cl68: return os << "68% C.L.";
      case Grid::header_t::cl95: return os << "95% C.L.";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const Grid::header_t::nucleon_t& nucl )
  {
    switch ( nucl ) {
      case Grid::header_t::proton: return os << "proton";
      case Grid::header_t::neutron: return os << "neutron";
    }
    return os;
  }
}

REGISTER_STRFUN( MSTWgrid, mstw::Grid )
