#ifndef CepGen_IO_MSTWGridHandler_h
#define CepGen_IO_MSTWGridHandler_h

#include <gsl/gsl_version.h>

#ifdef GSL_MAJOR_VERSION
#if GSL_MAJOR_VERSION > 2 || ( GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 1 )
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#define GOOD_GSL 1
#endif
#endif

#include <array>
#include <vector>
#include <set>

namespace CepGen { class StructureFunctions; }
namespace MSTW
{
  class GridHandler
  {
    public:
      struct sfval_t {
        float q2, xbj;
        double f2, fl;
      };
      struct header_t {
        enum order_t : unsigned short { lo = 0, nlo = 1, nnlo = 2 };
        enum cl_t : unsigned short { cl68 = 0, cl95 = 1 };
        enum nucleon_t : unsigned short { proton = 1, neutron = 2 };
        unsigned int magic;
        order_t order;
        cl_t cl;
        nucleon_t nucleon;
      };

    public:
      static GridHandler& get( const char* filename = "External/F2_Luxlike_fit/mstw_f2_scan_nnlo.dat" );
      ~GridHandler();

      CepGen::StructureFunctions eval( double q2, double xbj ) const;

      header_t header() const { return header_; }
      std::vector<sfval_t> values() const { return values_raw_; }

    private:
      explicit GridHandler( const char* );
      void initGSL( const std::set<double>& q2_vals, const std::set<double>& xbj_vals );

      enum spline_type { F2 = 0, FL = 1, num_functions_ };
      static const unsigned int good_magic;

      header_t header_;
      std::vector<sfval_t> values_raw_;
#ifdef GOOD_GSL
      std::array<gsl_spline2d*,2> splines_;
      gsl_interp_accel* xacc_, *yacc_;
      std::array<double*,2> values_;
#else
      std::vector<double> xbj_vals_, q2_vals_;
#endif

    public:
      GridHandler( const GridHandler& ) = delete;
      void operator=( const GridHandler& ) = delete;
  };
  std::ostream& operator<<( std::ostream&, const GridHandler::sfval_t& );
  std::ostream& operator<<( std::ostream&, const GridHandler::header_t::order_t& );
  std::ostream& operator<<( std::ostream&, const GridHandler::header_t::cl_t& );
  std::ostream& operator<<( std::ostream&, const GridHandler::header_t::nucleon_t& );
}

#endif

