#include "CepGen/Physics/GluonGrid.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/Timer.h"

#include <fstream>
#include <set>
#include <cmath>

namespace kmr
{
  GluonGrid&
  GluonGrid::get( const char* filename )
  {
    static GluonGrid instance( cepgen::ParametersList()
      .set<std::string>( "path", filename ) );
    return instance;
  }

  GluonGrid::GluonGrid( const cepgen::ParametersList& params ) :
    cepgen::GridHandler<3,1>( cepgen::GridType::linear ), // grid is already logarithmic
    grid_path_( params.get<std::string>( "path", DEFAULT_KMR_GRID_PATH ) )
  {
    CG_INFO( "GluonGrid" ) << "Building the KMR grid evaluator.";

    cepgen::utils::Timer tmr;

    std::set<double> kt2_vals, x_vals, mu2_vals;
    { // file readout part
      std::ifstream file( grid_path_, std::ios::in );
      if ( !file.is_open() )
        throw CG_FATAL( "GluonGrid" ) << "Failed to load grid file \"" << grid_path_ << "\"!";

      std::string line, x_tmp, kt2_tmp, mu2_tmp, fg_tmp;
      while ( file >> x_tmp >> kt2_tmp >> mu2_tmp >> fg_tmp ) {
        const double x = stod( x_tmp ), kt2 = stod( kt2_tmp ), mu2 = stod( mu2_tmp ), fg = stod( fg_tmp );
        x_vals.insert( x );
        kt2_vals.insert( kt2 );
        mu2_vals.insert( mu2 );
        insert( { x, kt2, mu2 }, { fg } );
      }
      file.close();
    }

    init();

    CG_INFO( "GluonGrid" )
      << "KMR grid evaluator built in " << tmr.elapsed() << " s.\n\t"
      << " kt² in range [" << *kt2_vals.begin() << ":" << *kt2_vals.rbegin() << "]\n\t"
      << " x in range [" << *x_vals.begin() << ":" << *x_vals.rbegin() << "]\n\t"
      << " µ² in range ["  << *mu2_vals.begin() << ":" << *mu2_vals.rbegin() << "].";
  }

  double
  GluonGrid::operator()( double x, double kt2, double mu2 ) const
  {
    return cepgen::GridHandler<3,1>::eval( { log( x ), log( kt2 ), log( mu2 ) } ).at( 0 );
  }
}
