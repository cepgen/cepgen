#ifndef CepGen_IO_MSTWGridHandler_h
#define CepGen_IO_MSTWGridHandler_h

#include "CepGen/IO/GridHandler.h"

namespace MSTW
{
  /// \note x/y = Q2/xbj given by the parent
  typedef CepGen::GridHandler<2>::grid_t sfval_t;
  std::ostream& operator<<( std::ostream&, const sfval_t& );

  class Grid : public CepGen::GridHandler<2>
  {
    public:
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
      static Grid& get( const char* filename = "External/F2_Luxlike_fit/mstw_f2_scan_nnlo.dat" );
      CepGen::StructureFunctions eval( double q2, double xbj ) const;
      header_t header() const { return header_; }

    public:
      Grid( const Grid& ) = delete;
      void operator=( const GridHandler& ) = delete;

    private:
      explicit Grid( const char* );

      enum spline_type { F2 = 0, FL = 1, num_functions_ };
      static const unsigned int good_magic;

      header_t header_;
  };
  std::ostream& operator<<( std::ostream&, const Grid::header_t::order_t& );
  std::ostream& operator<<( std::ostream&, const Grid::header_t::cl_t& );
  std::ostream& operator<<( std::ostream&, const Grid::header_t::nucleon_t& );
}

#endif

