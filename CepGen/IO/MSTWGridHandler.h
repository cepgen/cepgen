#ifndef CepGen_IO_MSTWGridHandler_h
#define CepGen_IO_MSTWGridHandler_h

#include <gsl/gsl_version.h>

#if GSL_MAJOR_VERSION > 2 || ( GSL_MAJOR_VERSION == 2 && GSL_MINOR_VERSION >= 1 )
#define GOOD_GSL 1
#endif

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#ifdef GOOD_GSL
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#endif

#include <fstream>
#include <array>
#include <vector>
#include <set>

namespace MSTW
{
  class GridHandler
  {
    public:
      static GridHandler& get( const char* filename );
      ~GridHandler();

      CepGen::StructureFunctions eval( double q2, double xbj ) const;

    private:
      GridHandler( const char* );

      struct sfval_t {
        float q2, xbj;
        double f2, fl;
      };
      struct header_t {
        unsigned int magic;
        enum { lo = 0, nlo = 1, nnlo = 2 } order;
      };
      static constexpr unsigned int good_magic = 0x5754534d; // MSTW in ASCII

#ifdef GOOD_GSL
      std::array<gsl_spline2d*,2> splines_;
      gsl_interp_accel* xacc_, *yacc_;
      std::array<double*,2> values_;
#endif

    public:
      GridHandler( const GridHandler& ) = delete;
      void operator=( const GridHandler& ) = delete;
  };
}

#endif

