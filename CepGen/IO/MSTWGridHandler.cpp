#include "MSTWGridHandler.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

namespace MSTW
{
  GridHandler&
  GridHandler::get( const char* filename )
  {
    static GridHandler instance( filename );
    return instance;
  }

  GridHandler::GridHandler( const char* filename )
#ifdef GOOD_GSL
    : xacc_( 0 ), yacc_( 0 ), values_( { { 0, 0 } } )
#endif
  {
    std::set<double> q2_vals, xbj_vals;

    { // file readout part
      std::ifstream file( filename, std::ios::binary | std::ios::in );
      if ( !file.is_open() )
        FatalError( Form( "Impossible to load grid file \"%s\"!", filename ) );

      file.read( reinterpret_cast<char*>( &header_ ), sizeof( header_t ) );

      // first checks on the file header

      if ( header_.magic != good_magic )
        FatalError( Form( "Wrong magic number retrieved: %u, expecting %u!", header_.magic, good_magic ) );

      if ( header_.nucleon != header_t::proton )
        FatalError( "Only proton structure function grids can be retrieved for this purpose!" );

      // retrieve all points and evaluate grid boundaries

      sfval_t val;
      while ( file.read( reinterpret_cast<char*>( &val ), sizeof( sfval_t ) ) ) {
        q2_vals.insert( log10( val.q2 ) );
        xbj_vals.insert( log10( val.xbj ) );
#ifndef GOOD_GSL
        q2_vals_.emplace_back( val.q2 );
        xbj_vals_.emplace_back( val.xbj );
#endif
        values_raw_.emplace_back( val );
      }
      file.close();
    }

    if ( q2_vals.size() < 2 || xbj_vals.size() < 2 )
      FatalError( "Invalid grid retrieved!" );

    initGSL( q2_vals, xbj_vals );

    {
      std::ostringstream ss_cl, ss_ord, ss_nucl;
      ss_cl << header_.cl;
      ss_ord << header_.order;
      ss_nucl << header_.nucleon;
      Information( Form( "MSTW@%s grid evaluator built for %s structure functions (%s)\n\t"
                         " Q² in range [%.3e:%.3e]\n\t"
                         "xBj in range [%.3e:%.3e]",
                         ss_ord.str().c_str(), ss_nucl.str().c_str(), ss_cl.str().c_str(),
                         pow( 10.,  *q2_vals.begin() ), pow( 10.,  *q2_vals.rbegin() ),
                         pow( 10., *xbj_vals.begin() ), pow( 10., *xbj_vals.rbegin() ) ) );
    }
  }

  GridHandler::~GridHandler()
  {
#ifdef GOOD_GSL
    for ( unsigned short i = 0; i < num_functions_; ++i ) {
      gsl_spline2d_free( splines_[i] );
      if ( values_[i] ) delete[] values_[i];
    }
    if ( xacc_ ) gsl_interp_accel_free( xacc_ );
    if ( yacc_ ) gsl_interp_accel_free( yacc_ );
#endif
  }

  void
  GridHandler::initGSL( const std::set<double>& q2_vals, const std::set<double>& xbj_vals )
  {
#ifdef GOOD_GSL
    gsl_set_error_handler_off();
    const gsl_interp2d_type* T = gsl_interp2d_bilinear;

    for ( unsigned short i = 0; i < num_functions_; ++i ) {
      values_[i] = new double[q2_vals.size() * xbj_vals.size()];
      splines_[i] = gsl_spline2d_alloc( T, q2_vals.size(), xbj_vals.size() );
    }
    xacc_ = gsl_interp_accel_alloc();
    yacc_ = gsl_interp_accel_alloc();

    // second loop over all points to populate the grid
    for ( const auto& val : values_raw_ ) {
      // retrieve the index of the Q2/xbj bin in the set
      const unsigned short id_q2  = std::distance(  q2_vals.begin(),  q2_vals.lower_bound( log10( val.q2  ) ) ),
                           id_xbj = std::distance( xbj_vals.begin(), xbj_vals.lower_bound( log10( val.xbj ) ) );
      gsl_spline2d_set( splines_[F2], values_[F2], id_q2, id_xbj, val.f2 );
      gsl_spline2d_set( splines_[FL], values_[FL], id_q2, id_xbj, val.fl );
    }

    // initialise splines objects
    std::vector<double> q2_vec( q2_vals.begin(), q2_vals.end() ), xbj_vec( xbj_vals.begin(), xbj_vals.end() );
    for ( unsigned short i = 0; i < num_functions_; ++i ) {
      gsl_spline2d_init( splines_[i], &q2_vec[0], &xbj_vec[0], values_[i], q2_vals.size(), xbj_vals.size() );
    }
#else
    InWarning( Form( "GSL version ≥ 2.1 is required for spline bilinear interpolation.\n\t"
                     "Version %s is installed on this system!\n\t"
                     "Will use a linear approximation instead.\n\t"
                     "You may check the numerical validity of this approach...", GSL_VERSION ) );
#endif
  }

  CepGen::StructureFunctions
  GridHandler::eval( double q2, double xbj ) const
  {
    CepGen::StructureFunctions ev;
#ifdef GOOD_GSL
    if ( gsl_spline2d_eval_e( splines_[F2], log10( q2 ), log10( xbj ), xacc_, yacc_, &ev.F2 ) != GSL_SUCCESS
      || gsl_spline2d_eval_e( splines_[FL], log10( q2 ), log10( xbj ), xacc_, yacc_, &ev.FL ) != GSL_SUCCESS ) {
      InWarning( Form( "Failed to evaluate the structure functions for Q² = %.5e GeV² / xbj = %.5e", q2, xbj ) );
      return ev;
    }
#else
    double q2_1 = -1., q2_2 = -1., xbj_1 = -1., xbj_2 = -1.;
    double f2_11 = 0., f2_12 = 0., f2_21 = 0., f2_22 = 0.;
    double fl_11 = 0., fl_12 = 0., fl_21 = 0., fl_22 = 0.;
    for ( unsigned int i = 0; i < q2_vals_.size()-1; ++i ) {
      if ( q2 >= q2_vals_.at( i ) && q2 < q2_vals_.at( i+1 ) ) {
        q2_1 = q2_vals_.at( i );
        q2_2 = q2_vals_.at( i+1 );
        break;
      }
    }
    for ( unsigned int i = 0; i < xbj_vals_.size()-1; ++i ) {
      if ( xbj >= xbj_vals_.at( i ) && xbj < xbj_vals_.at( i+1 ) ) {
        xbj_1 = xbj_vals_.at( i );
        xbj_2 = xbj_vals_.at( i+1 );
        break;
      }
    }
    for ( const auto& val : values_raw_ ) {
      if ( q2_1 == val.q2 && xbj_1 == val.xbj ) {
        f2_11 = val.f2;
        fl_11 = val.fl;
      }
      if ( q2_2 == val.q2 && xbj_1 == val.xbj ) {
        f2_12 = val.f2;
        fl_12 = val.fl;
      }
      if ( q2_1 == val.q2 && xbj_2 == val.xbj ) {
        f2_21 = val.f2;
        fl_21 = val.fl;
      }
      if ( q2_2 == val.q2 && xbj_2 == val.xbj ) {
        f2_22 = val.f2;
        fl_22 = val.fl;
      }
    }
    const double x2x1 = xbj_2-xbj_1;
    const double y2y1 = q2_2-q2_1;
    const double x2x = xbj_2-xbj;
    const double y2y = q2_2-q2;
    const double yy1 = q2-q2_1;
    const double xx1 = xbj-xbj_1;
    ev.F2 = 1.0 / (x2x1 * y2y1) * ( f2_11*x2x*y2y + f2_21*xx1*y2y + f2_12*x2x*yy1 + f2_22*xx1*yy1 );
    ev.FL = 1.0 / (x2x1 * y2y1) * ( fl_11*x2x*y2y + fl_21*xx1*y2y + fl_12*x2x*yy1 + fl_22*xx1*yy1 );
#endif
    return ev;
  }

  std::ostream&
  operator<<( std::ostream& os, const GridHandler::sfval_t& val )
  {
    return os << Form( "Q² = %.5e GeV²\txbj = %.4f\tF₂ = % .6e\tFL = % .6e", val.q2, val.xbj, val.f2, val.fl );
  }

  std::ostream&
  operator<<( std::ostream& os, const GridHandler::header_t::order_t& order )
  {
    switch ( order ) {
      case GridHandler::header_t::lo: return os << "LO";
      case GridHandler::header_t::nlo: return os << "nLO";
      case GridHandler::header_t::nnlo: return os << "nnLO";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const GridHandler::header_t::cl_t& cl )
  {
    switch ( cl ) {
      case GridHandler::header_t::cl68: return os << "68% C.L.";
      case GridHandler::header_t::cl95: return os << "95% C.L.";
    }
    return os;
  }

  std::ostream&
  operator<<( std::ostream& os, const GridHandler::header_t::nucleon_t& nucl )
  {
    switch ( nucl ) {
      case GridHandler::header_t::proton: return os << "proton";
      case GridHandler::header_t::neutron: return os << "neutron";
    }
    return os;
  }
}
