/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2025  Laurent Forthomme
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <fstream>

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/GridHandler.h"
#include "CepGen/Utils/String.h"

/// Martin-Stirling-Thorne-Watt PDFs structure functions
namespace mstw {
  /// A \f$F_{2,L}\f$ grid interpolator
  class Grid final : public cepgen::strfun::Parameterisation, cepgen::GridHandler<2, 2> {
  public:
    /// Grid MSTW structure functions evaluator
    explicit Grid(const cepgen::ParametersList& params)
        : Parameterisation(params), GridHandler(cepgen::GridType::logarithmic) {
      {  // file readout part
        const auto& grid_path = steerPath("gridPath");
        std::ifstream file(grid_path, std::ios::binary | std::ios::in);
        if (!file.is_open())
          throw CG_FATAL("MSTW") << "Failed to load grid file \"" << grid_path << "\"!";

        file.read(reinterpret_cast<char*>(&header_), sizeof(header_t));

        // first checks on the file header

        if (header_.magic != GOOD_MAGIC)
          throw CG_FATAL("MSTW") << "Wrong magic number retrieved: " << header_.magic << ", expecting " << GOOD_MAGIC
                                 << ".";

        if (header_.nucleon != header_t::proton)
          throw CG_FATAL("MSTW") << "Only proton structure function grids can be retrieved for this purpose!";

        // retrieve all points and evaluate grid boundaries

        sfval_t val{};
        while (file.read(reinterpret_cast<char*>(&val), sizeof(sfval_t)))
          insert({val.xbj, val.q2}, {val.f2, val.fl});
        file.close();
        initialise();  // initialise the grid after filling its nodes
      }
      const auto& bounds = boundaries();
      CG_DEBUG("MSTW") << "MSTW@" << header_.order << " grid evaluator built "
                       << "for " << header_.nucleon << " structure functions (" << header_.cl << ")\n\t"
                       << "xBj in range [" << std::pow(10., bounds[0].min()) << ":" << std::pow(10., bounds[0].max())
                       << "], Q² in range [" << std::pow(10., bounds[1].min()) << ":" << std::pow(10., bounds[1].max())
                       << "].";
    }

    static cepgen::ParametersDescription description() {
      auto desc = Parameterisation::description();
      desc.setDescription("MSTW grid (perturbative)");
      desc.add<std::string>("gridPath", "mstw_sf_scan_nnlo.dat").setDescription("Path to the MSTW grid content");
      return desc;
    }

    /// Grid header information as parsed from the file
    struct header_t {
      /// Interpolation order
      enum order_t : unsigned short { lo = 0, nlo = 1, nnlo = 2 };
      /// Confidence level
      enum cl_t : unsigned short { cl68 = 0, cl95 = 1 };
      /// Type of nucleon interpolated
      enum nucleon_t : unsigned short { proton = 1, neutron = 2 };
      unsigned int magic;  ///< Grid file magic number
      order_t order;       ///< Interpolation order
      cl_t cl;             ///< Confidence level
      nucleon_t nucleon;   ///< Type of nucleon interpolated
    };
    /// Structure functions value at a given \f$Q^2/x_{\rm Bj}\f$ coordinate
    struct sfval_t {
      float q2{0.};   ///< four-momentum transfer, in GeV\f$^2\f$
      float xbj{0.};  ///< Bjorken's scaling variable
      double f2{0.};  ///< Transverse structure function value
      double fl{0.};  ///< Longitudinal structure function value
    };

    /// Compute the structure functions at a given \f$Q^2/x_{\rm Bj}\f$
    void eval() override {
      const auto& val = cepgen::GridHandler<2, 2>::eval({args_.xbj, args_.q2});
      setF2(val.at(0));
      setFL(val.at(1));
    }
    /// Retrieve the grid's header information
    header_t header() const { return header_; }

    //--- already retrieved from grid, so no need to recompute it
    Grid& computeFL(double, double) override { return *this; }
    Grid& computeFL(double, double, double) override { return *this; }

  private:
    static constexpr unsigned int GOOD_MAGIC = 0x5754534d;  // MSTW in ASCII

    header_t header_ = {};
  };

  /// Human-readable description of a values point
  std::ostream& operator<<(std::ostream& os, const Grid::sfval_t& val) {
    return os << cepgen::utils::format(
               "xbj = %.4f\tQ² = %.5e GeV²\tF_2 = % .6e\tF_1 = % .6e", val.xbj, val.q2, val.f2, val.fl);
  }

  /// Human-readable description of an interpolation order
  std::ostream& operator<<(std::ostream& os, const Grid::header_t::order_t& order) {
    switch (order) {
      case Grid::header_t::lo:
        return os << "LO";
      case Grid::header_t::nlo:
        return os << "nLO";
      case Grid::header_t::nnlo:
        return os << "nnLO";
    }
    return os;
  }

  /// Human-readable description of a confidence level
  std::ostream& operator<<(std::ostream& os, const Grid::header_t::cl_t& cl) {
    switch (cl) {
      case Grid::header_t::cl68:
        return os << "68% C.L.";
      case Grid::header_t::cl95:
        return os << "95% C.L.";
    }
    return os;
  }

  /// Human-readable description of a nucleon type
  std::ostream& operator<<(std::ostream& os, const Grid::header_t::nucleon_t& nucleon) {
    switch (nucleon) {
      case Grid::header_t::proton:
        return os << "proton";
      case Grid::header_t::neutron:
        return os << "neutron";
    }
    return os;
  }
}  // namespace mstw
using mstw::Grid;
REGISTER_STRFUN("MSTWGrid", 205, Grid);
