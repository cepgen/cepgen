#ifndef CepGen_StructureFunctions_MSTWGrid_h
#define CepGen_StructureFunctions_MSTWGrid_h

#include "CepGen/IO/GridHandler.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#define DEFAULT_MSTW_GRID_PATH "External/mstw_sf_scan_nnlo.dat"

/// Martin-Stirling-Thorne-Watt PDFs structure functions
namespace MSTW
{
  /// A \f$F_{2,L}\f$ grid interpolator
  class Grid : public CepGen::StructureFunctions, private CepGen::GridHandler<2,2>
  {
    public:
      /// Grid header information as parsed from the file
      struct header_t
      {
        enum order_t : unsigned short { lo = 0, nlo = 1, nnlo = 2 };
        enum cl_t : unsigned short { cl68 = 0, cl95 = 1 };
        enum nucleon_t : unsigned short { proton = 1, neutron = 2 };
        unsigned int magic;
        order_t order;
        cl_t cl;
        nucleon_t nucleon;
      };
      struct sfval_t
      {
        float q2, xbj;
        double f2, fl;
      };
      struct Parameterisation {
        Parameterisation() : grid_path( DEFAULT_MSTW_GRID_PATH ) {}
        std::string grid_path;
      };

    public:
      /// Retrieve the grid interpolator (singleton)
      static Grid& get( const char* path = DEFAULT_MSTW_GRID_PATH );

      /// Compute the structure functions at a given \f$Q^2/x_{\rm Bj}\f$
      Grid& operator()( double xbj, double q2 ) override;
      /// Retrieve the grid's header information
      header_t header() const { return header_; }
      Parameterisation params;

        //--- already retrieved from grid, so no need to recompute it
      void computeFL( double xbj, double q2, const CepGen::SF::SigmaRatio& ) override {}
      void computeFL( double xbj, double q2, double r ) override {}

    public:
      Grid( const Grid& ) = delete;
      void operator=( const GridHandler& ) = delete;

    private:
      explicit Grid( const Parameterisation& = Parameterisation() );
      static const unsigned int good_magic;
      static std::shared_ptr<Grid> singl_;

      header_t header_;
  };

  std::ostream& operator<<( std::ostream&, const Grid::sfval_t& );
  std::ostream& operator<<( std::ostream&, const Grid::header_t::order_t& );
  std::ostream& operator<<( std::ostream&, const Grid::header_t::cl_t& );
  std::ostream& operator<<( std::ostream&, const Grid::header_t::nucleon_t& );
}

#undef DEFAULT_MSTW_GRID_PATH

#endif

