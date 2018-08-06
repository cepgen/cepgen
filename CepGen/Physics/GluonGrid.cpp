#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Core/Exception.h"

#include <fstream>
#include <set>

namespace kmr
{
  GluonGrid&
  GluonGrid::get( const char* filename )
  {
    Parameterisation p;
    p.grid_path = filename;
    static GluonGrid instance( p );
    return instance;
  }

  GluonGrid::GluonGrid( const Parameterisation& param ) :
    CepGen::GridHandler<3,1>( CepGen::GridType::linear ),
    params( param )
  {
    std::set<double> kt2_vals, x_vals, mu2_vals;

    { // file readout part
      std::ifstream file( params.grid_path, std::ios::in );
      if ( !file.is_open() )
        throw CG_FATAL( "GluonGrid" ) << "Impossible to load grid file \"" << params.grid_path << "\"!";

      std::string x_tmp, kt2_tmp, mu2_tmp, fg_tmp;
      while ( file >> x_tmp >> kt2_tmp >> mu2_tmp >> fg_tmp ) {
        const double x = stod( x_tmp ), kt2 = stod( kt2_tmp ), mu2 = stod( mu2_tmp ), fg = stod( fg_tmp );
        kt2_vals.insert( kt2 );
        x_vals.insert( x );
        mu2_vals.insert( mu2 );
        insert( CepGen::GridHandler<3,1>::point_t{ { kt2, x, mu2 }, { fg } } );
      }
      file.close();
    }

    init();

    CG_INFO( "GluonGrid" )
      << "KMR grid evaluator built!\n\t"
      << " kt² in range [" << *kt2_vals.begin() << ":" << *kt2_vals.rbegin() << "]\n\t"
      << " x in range [" << *x_vals.begin() << ":" << *x_vals.rbegin() << "]\n\t"
      << " µ² in range ["  << *mu2_vals.begin() << ":" << *mu2_vals.rbegin() << "].";
  }

  double
  GluonGrid::operator()( double kt2, double x, double mu2 ) const
  {
    return CepGen::GridHandler<3,1>::eval( { kt2, x, mu2 } ).at( 0 );
  }
}

