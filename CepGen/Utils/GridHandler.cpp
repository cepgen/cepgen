#include "CepGen/Utils/GridHandler.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <limits>

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

namespace cepgen {
  template <size_t D, size_t N>
  GridHandler<D, N>::GridHandler(const GridType& grid_type) : grid_type_(grid_type), accel_{}, init_(false) {
    for (size_t i = 0; i < D; ++i)
      accel_.emplace_back(gsl_interp_accel_alloc(), gsl_interp_accel_free);
  }

  template <size_t D, size_t N>
  typename GridHandler<D, N>::values_t GridHandler<D, N>::eval(coord_t in_coords) const {
    if (!init_)
      throw CG_FATAL("GridHandler") << "Grid extrapolator called but not initialised!";

    values_t out;
    coord_t coord = in_coords;
    switch (grid_type_) {
      case GridType::logarithmic: {
        for (auto& c : coord)
          c = log10(c);
      } break;
      case GridType::square: {
        for (auto& c : coord)
          c *= c;
      } break;
      default:
        break;
    }
    //--- dimension of the vector space coordinate to evaluate
    switch (D) {
      case 1: {
        for (size_t i = 0; i < N; ++i) {
          int res = gsl_spline_eval_e(splines_1d_.at(i).get(), coord.at(0), accel_.at(0).get(), &out[i]);
          if (res != GSL_SUCCESS) {
            out[i] = 0.;
            CG_WARNING("GridHandler") << "Failed to evaluate the grid value (N=" << i << ") "
                                      << "for x = " << in_coords.at(0) << ". "
                                      << "GSL error: " << gsl_strerror(res);
          }
        }
      } break;
      case 2: {
#ifdef GSL_VERSION_ABOVE_2_1
        const double x = coord.at(0), y = coord.at(1);
        for (size_t i = 0; i < N; ++i) {
          int res = gsl_spline2d_eval_e(splines_2d_.at(i).get(), x, y, accel_.at(0).get(), accel_.at(1).get(), &out[i]);
          if (res != GSL_SUCCESS) {
            out[i] = 0.;
            CG_WARNING("GridHandler") << "Failed to evaluate the grid value (N=" << i << ") "
                                      << "for x = " << x << " / y = " << y << ". "
                                      << "GSL error: " << gsl_strerror(res);
          }
        }
#else
        //--- retrieve the indices of the bin in the set
        coord_t before(D), after(D);
        findIndices(coord, before, after);
        //--- find boundaries values
        const gridpoint_t &ext_11 = values_raw_.at({before[0], before[1]}),
                          &ext_12 = values_raw_.at({before[0], after[1]}),
                          &ext_21 = values_raw_.at({after[0], before[1]}),
                          &ext_22 = values_raw_.at({after[0], after[1]});
        //--- now that we have the boundaries, we may interpolate
        coord_t c_d(D);
        for (size_t i = 0; i < D; ++i)
          c_d[i] = (after[i] != before[i]) ? (coord.at(i) - before[i]) / (after[i] - before[i]) : 0.;
        const gridpoint_t ext_1 = ext_11 * (1. - c_d[0]) + ext_21 * c_d[0];
        const gridpoint_t ext_2 = ext_12 * (1. - c_d[0]) + ext_22 * c_d[0];
        out = ext_1 * (1. - c_d[1]) + ext_2 * c_d[1];
#endif
      } break;
      case 3: {
        //--- retrieve the indices of the bin in the set
        coord_t before(D), after(D);
        findIndices(coord, before, after);
        //--- find boundaries values
        const gridpoint_t &ext_111 = values_raw_.at({before[0], before[1], before[2]}),
                          &ext_112 = values_raw_.at({before[0], before[1], after[2]}),
                          &ext_121 = values_raw_.at({before[0], after[1], before[2]}),
                          &ext_122 = values_raw_.at({before[0], after[1], after[2]}),
                          &ext_211 = values_raw_.at({after[0], before[1], before[2]}),
                          &ext_212 = values_raw_.at({after[0], before[1], after[2]}),
                          &ext_221 = values_raw_.at({after[0], after[1], before[2]}),
                          &ext_222 = values_raw_.at({after[0], after[1], after[2]});
        //--- now that we have the boundaries, we may interpolate
        coord_t c_d(D);
        for (size_t i = 0; i < D; ++i)
          c_d[i] = (after[i] != before[i]) ? (coord.at(i) - before[i]) / (after[i] - before[i]) : 0.;
        const gridpoint_t ext_11 = ext_111 * (1. - c_d[0]) + ext_211 * c_d[0];
        const gridpoint_t ext_12 = ext_112 * (1. - c_d[0]) + ext_212 * c_d[0];
        const gridpoint_t ext_21 = ext_121 * (1. - c_d[0]) + ext_221 * c_d[0];
        const gridpoint_t ext_22 = ext_122 * (1. - c_d[0]) + ext_222 * c_d[0];
        const gridpoint_t ext_1 = ext_11 * (1. - c_d[1]) + ext_21 * c_d[1];
        const gridpoint_t ext_2 = ext_12 * (1. - c_d[1]) + ext_22 * c_d[1];
        out = ext_1 * (1. - c_d[2]) + ext_2 * c_d[2];
      } break;
      default:
        throw CG_FATAL("GridHandler") << "Unsupported number of dimensions: " << N << ".\n\t"
                                      << "Please contact the developers to add such a new feature.";
    }
    return out;
  }

  template <size_t D, size_t N>
  void GridHandler<D, N>::insert(coord_t coord, values_t value) {
    auto mod_coord = coord;
    if (grid_type_ != GridType::linear)
      for (auto& c : mod_coord)
        switch (grid_type_) {
          case GridType::logarithmic:
            c = log10(c);
            break;
          case GridType::square:
            c *= c;
            break;
          default:
            break;
        }
    if (values_raw_.count(mod_coord) != 0)
      CG_WARNING("GridHandler") << "Duplicate coordinate detected for x=" << coord << ".";
    values_raw_[mod_coord] = value;
    init_ = false;
  }

  template <size_t D, size_t N>
  void GridHandler<D, N>::init() {
    if (values_raw_.empty())
      throw CG_ERROR("GridHandler") << "Empty grid.";
    gsl_set_error_handler_off();
    //--- start by building grid coordinates from raw values
    for (auto& c : coords_)
      c.clear();
    for (const auto& val : values_raw_) {
      unsigned short i = 0;
      for (const auto& c : val.first) {
        if (std::find(coords_.at(i).begin(), coords_.at(i).end(), c) == coords_.at(i).end())
          coords_.at(i).emplace_back(c);
        ++i;
      }
    }
    for (auto& c : coords_)
      std::sort(c.begin(), c.end());
    CG_DEBUG("GridHandler").log([&](auto& dbg) {
      dbg << "Grid dump:";
      // debugging of the grid coordinates
      unsigned short i = 0;
      for (const auto& cs : coords_) {
        dbg << "\n>> coordinate " << (i++) << " has " << utils::s("member", cs.size(), true) << ":";
        unsigned short j = 0;
        for (const auto& val : cs)
          dbg << (j++ % 20 == 0 ? "\n  " : " ") << val;
      }
    });
    //--- particularise by dimension
    switch (D) {
      case 1: {  //--- x |-> (f1,...)
        const gsl_interp_type* type = gsl_interp_cspline;
        //const gsl_interp_type* type = gsl_interp_steffen;
#ifdef GSL_VERSION_ABOVE_2_1
        const unsigned short min_size = gsl_interp_type_min_size(type);
#else
        const unsigned short min_size = type->min_size;
#endif
        if (min_size >= values_raw_.size())
          throw CG_FATAL("GridHandler") << "Not enough points for \"" << type->name << "\" type of interpolation.\n\t"
                                        << "Minimum required: " << min_size << ", got " << values_raw_.size() << "!";
        for (size_t i = 0; i < N; ++i) {
          values_[i].reset(new double[values_raw_.size()]);
          splines_1d_.emplace_back(gsl_spline_alloc(type, values_raw_.size()), gsl_spline_free);
        }
        // transform a map(coord -> values) to
        // - 1 vector(coordinates)
        // - N vectors(values)
        std::vector<double> x_vec;
        size_t i = 0;
        for (const auto& vals : values_raw_) {
          x_vec.emplace_back(vals.first.at(0));
          unsigned short j = 0;
          for (const auto& val : vals.second)
            values_[j++].get()[i] = val;
          ++i;
        }
        // initialise spline interpolation objects (one for each value)
        for (size_t j = 0; j < splines_1d_.size(); ++j)
          gsl_spline_init(splines_1d_.at(j).get(), &x_vec[0], values_[j].get(), values_raw_.size());
      } break;
      case 2: {  //--- (x,y) |-> (f1,...)
#ifdef GSL_VERSION_ABOVE_2_1
        const gsl_interp2d_type* type = gsl_interp2d_bilinear;
        splines_2d_.clear();
        for (size_t i = 0; i < N; ++i) {
          values_[i].reset(new double[coords_.at(0).size() * coords_.at(1).size()]);
          splines_2d_.emplace_back(gsl_spline2d_alloc(type, coords_.at(0).size(), coords_.at(1).size()),
                                   gsl_spline2d_free);
        }

        // second loop over all points to populate the grid
        for (const auto& val : values_raw_) {
          double val_x = val.first.at(0), val_y = val.first.at(1);
          // retrieve the index of the bin in the set
          const unsigned short id_x =
              std::distance(coords_.at(0).begin(), std::lower_bound(coords_.at(0).begin(), coords_.at(0).end(), val_x));
          const unsigned short id_y =
              std::distance(coords_.at(1).begin(), std::lower_bound(coords_.at(1).begin(), coords_.at(1).end(), val_y));
          for (unsigned short i = 0; i < splines_2d_.size(); ++i)
            gsl_spline2d_set(splines_2d_.at(i).get(), values_[i].get(), id_x, id_y, val.second[i]);
        }

        // initialise spline interpolation objects (one for each value)
        const coord_t &x_vec = coords_.at(0), &y_vec = coords_.at(1);
        for (size_t i = 0; i < splines_2d_.size(); ++i)
          gsl_spline2d_init(
              splines_2d_.at(i).get(), &x_vec[0], &y_vec[0], values_[i].get(), x_vec.size(), y_vec.size());
#else
        CG_WARNING("GridHandler") << "GSL version ≥ 2.1 is required for spline bilinear interpolation.\n\t"
                                  << "Version " << GSL_VERSION << " is installed on this system!\n\t"
                                  << "Will use a simple bilinear approximation instead.";
#endif
      } break;
      default:
        break;
    }
    init_ = true;
    CG_DEBUG("GridHandler") << "Grid evaluator initialised with boundaries: " << boundaries() << "\n"
                            << "Values handled:\n"
                            << values_raw_;
  }

  template <size_t D, size_t N>
  std::array<std::pair<double, double>, D> GridHandler<D, N>::boundaries() const {
    std::array<std::pair<double, double>, D> out;
    const auto min_val = min(), max_val = max();
    for (size_t i = 0; i < D; ++i)
      out[i] = {min_val[i], max_val[i]};
    return out;
  }

  template <size_t D, size_t N>
  std::array<double, D> GridHandler<D, N>::min() const {
    std::array<double, D> out;
    size_t i = 0;
    for (const auto& c : coords_) {  // loop over all dimensions
      const auto& min = std::min_element(c.begin(), c.end());
      out[i++] = min != c.end() ? *min : std::numeric_limits<double>::infinity();
    }
    return out;
  }

  template <size_t D, size_t N>
  std::array<double, D> GridHandler<D, N>::max() const {
    std::array<double, D> out;
    size_t i = 0;
    for (const auto& c : coords_) {  // loop over all dimensions
      const auto& max = std::max_element(c.begin(), c.end());
      out[i++] = max != c.end() ? *max : std::numeric_limits<double>::infinity();
    }
    return out;
  }

  template <size_t D, size_t N>
  void GridHandler<D, N>::findIndices(const coord_t& coord, coord_t& min, coord_t& max) const {
    for (size_t i = 0; i < D; ++i) {
      const auto& c_i = coords_.at(i);   // extract all coordinates registered for this dimension
      if (coord.at(i) < *c_i.begin()) {  // under the range
        CG_DEBUG_LOOP("GridHandler:indices") << "Coordinate " << i << " in underflow range "
                                             << "(" << coord.at(i) << " < " << *c_i.begin() << ").";
        min[i] = max[i] = *c_i.begin();
      } else if (coord.at(i) > *c_i.rbegin()) {  // over the range
        CG_DEBUG_LOOP("GridHandler:indices") << "Coordinate " << i << " in overflow range "
                                             << "(" << coord.at(i) << " > " << *c_i.rbegin() << ").";
        min[i] = max[i] = *c_i.rbegin();
      } else {  // in between two coordinates
        auto it_coord = std::lower_bound(c_i.begin(), c_i.end(), coord.at(i));
        max[i] = *it_coord;
        if (it_coord != c_i.begin())
          it_coord--;
        min[i] = *it_coord;
      }
    }
  }

  //----------------------------------------------------------------------------
  // grid manipulation utility
  //----------------------------------------------------------------------------

  template <size_t D, size_t N>
  typename GridHandler<D, N>::gridpoint_t GridHandler<D, N>::gridpoint_t::operator*(double c) const {
    gridpoint_t out = *this;
    for (auto& a : out)
      a *= c;
    return out;
  }

  template <size_t D, size_t N>
  typename GridHandler<D, N>::gridpoint_t GridHandler<D, N>::gridpoint_t::operator+(const gridpoint_t& rhs) const {
    gridpoint_t out = *this;
    for (size_t i = 0; i < out.size(); ++i)
      out[i] += rhs[i];
    return out;
  }

  template class GridHandler<1, 1>;
  template class GridHandler<1, 2>;
  template class GridHandler<2, 2>;
  template class GridHandler<3, 1>;
}  // namespace cepgen
