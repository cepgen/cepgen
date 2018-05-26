#include "CepGen/StructureFunctions/MSTWGrid.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

namespace MSTW
{
  const unsigned int Grid::good_magic = 0x5754534d; // MSTW in ASCII

  Grid&
  Grid::get( const char* filename )
  {
    static Grid instance( filename );
    return instance;
  }

  Grid::Grid( const char* filename )
  {
    std::set<double> q2_vals, xbj_vals;

    { // file readout part
      std::ifstream file( filename, std::ios::binary | std::ios::in );
      if ( !file.is_open() )
        throw CG_FATAL( "Grid" ) << "Impossible to load grid file \"" << filename << "%s\"!";

      file.read( reinterpret_cast<char*>( &header_ ), sizeof( header_t ) );

      // first checks on the file header

      if ( header_.magic != good_magic )
        throw CG_FATAL( "Grid" ) << "Wrong magic number retrieved: " << header_.magic << ", expecting " << good_magic << ".";

      if ( header_.nucleon != header_t::proton )
        throw CG_FATAL( "Grid" ) << "Only proton structure function grids can be retrieved for this purpose!";

      // retrieve all points and evaluate grid boundaries

      sfval_t val;
      while ( file.read( reinterpret_cast<char*>( &val ), sizeof( sfval_t ) ) ) {
        q2_vals.insert( log10( val.x ) );
        xbj_vals.insert( log10( val.y ) );
#ifndef GOOD_GSL
        x_vals_.emplace_back( val.x );
        y_vals_.emplace_back( val.y );
#endif
        values_raw_.emplace_back( val );
      }
      file.close();
    }

    if ( q2_vals.size() < 2 || xbj_vals.size() < 2 )
      throw CG_FATAL( "Grid" ) << "Invalid grid retrieved!";

    initGSL( q2_vals, xbj_vals );

    CG_INFO( "Grid" )
      << "MSTW@" << header_.order << " grid evaluator built "
      << "for " << header_.nucleon << " structure functions (" << header_.cl << ")\n\t"
      << " Q² in range [" << pow( 10., *q2_vals.begin() ) << ":" << pow( 10., *q2_vals.rbegin() ) << "]\n\t"
      << "xBj in range [" << pow( 10., *xbj_vals.begin() ) << ":" << pow( 10., *xbj_vals.rbegin() ) << "].";
  }

  CepGen::StructureFunctions
  Grid::eval( double q2, double xbj ) const
  {
    const sfval_t val = CepGen::GridHandler<2>::eval( q2, xbj );
    CepGen::StructureFunctions sf;
    sf.F2 = val.value[0];
    sf.FL = val.value[1];
    return sf;
  }

  std::ostream&
  operator<<( std::ostream& os, const sfval_t& val )
  {
    return os << Form( "Q² = %.5e GeV²\txbj = %.4f\tF₂ = % .6e\tFL = % .6e", val.x, val.y, val.value[0], val.value[1] );
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
