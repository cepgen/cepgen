#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <fstream>
#include <array>
#include <vector>
#include <set>

namespace MSTW
{
  class GridHandler
  {
    public:
      GridHandler( const char* filename );
      ~GridHandler();

      CepGen::StructureFunctions eval( double q2, double xbj ) const;

    private:
      std::array<gsl_spline2d*,2> splines_;
      gsl_interp_accel* xacc_, *yacc_;
      double* values_;
  };
}
