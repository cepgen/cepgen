#ifndef CepGen_StructureFunctions_MSTWGrid_h
#define CepGen_StructureFunctions_MSTWGrid_h

#include "CepGen/IO/GridHandler.h"
#include "CepGen/StructureFunctions/StructureFunctions.h"

#define DEFAULT_MSTW_GRID_PATH "External/mstw_sf_scan_nnlo.dat"

/// Martin-Stirling-Thorne-Watt PDFs structure functions
namespace mstw
{
  /// A \f$F_{2,L}\f$ grid interpolator
  class Grid : public CepGen::sf::Parameterisation, private CepGen::GridHandler<2,2>
  {
    public:
      /// Grid header information as parsed from the file
      struct header_t
      {
        /// Interpolation order
        enum order_t : unsigned short { lo = 0, nlo = 1, nnlo = 2 };
        /// Confidence level
        enum cl_t : unsigned short { cl68 = 0, cl95 = 1 };
        /// Type of nucleon interpolated
        enum nucleon_t : unsigned short { proton = 1, neutron = 2 };
        unsigned int magic; ///< Grid file magic number
        order_t order; ///< Interpolation order
        cl_t cl; ///< Confidence level
        nucleon_t nucleon; ///< Type of nucleon interpolated
      };
      /// Structure functions value at a given \f$Q^2/x_{\rm Bj}\f$ coordinate
      struct sfval_t
      {
        float q2; ///< four-momentum transfer, in GeV\f${}^2\f$
        float xbj; ///< Bjorken's scaling variable
        double f2; ///< Transverse structure function value
        double fl; ///< Longitudinal structure function value
      };
      /// List of parameters for this grid definition
      struct Parameters {
        Parameters() : grid_path( DEFAULT_MSTW_GRID_PATH ) {}
        /// Location of the grid to be interpolated
        std::string grid_path;
      };

    public:
      /// Retrieve the grid interpolator (singleton)
      static Grid& get( const char* path = DEFAULT_MSTW_GRID_PATH );

      /// Compute the structure functions at a given \f$Q^2/x_{\rm Bj}\f$
      Grid& operator()( double xbj, double q2 ) override;
      /// Retrieve the grid's header information
      header_t header() const { return header_; }
      /// Grid parameterisation object
      Parameters params;

        //--- already retrieved from grid, so no need to recompute it
      void computeFL( double xbj, double q2, const CepGen::sr::Parameterisation& ) override {}
      void computeFL( double xbj, double q2, double r ) override {}

    public:
      Grid( const Grid& ) = delete;
      void operator=( const GridHandler& ) = delete;

    private:
      explicit Grid( const Parameters& = Parameters() );
      std::string description() const override;
      static const unsigned int good_magic;
      static std::shared_ptr<Grid> singl_;

      header_t header_;
  };

  std::ostream& operator<<( std::ostream&, const Grid::sfval_t& ); ///< Human-readable description of a values point
  std::ostream& operator<<( std::ostream&, const Grid::header_t::order_t& ); ///< Human-readable description of an interpolation order
  std::ostream& operator<<( std::ostream&, const Grid::header_t::cl_t& ); ///< Human-readable description of a confidence level
  std::ostream& operator<<( std::ostream&, const Grid::header_t::nucleon_t& ); ///< Human-readable description of a nucleon type
}

#undef DEFAULT_MSTW_GRID_PATH

#endif
