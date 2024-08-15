/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2018-2024  Laurent Forthomme
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

#ifndef CepGen_Utils_Limits_h
#define CepGen_Utils_Limits_h

#include <iosfwd>
#include <utility>
#include <vector>

namespace cepgen {
  /// Validity interval for a variable
  class Limits : private std::pair<double, double> {
  public:
    Limits(double min = INVALID, double max = INVALID);  ///< Define lower and upper limits on a quantity
    Limits(const Limits&);                               ///< Copy constructor

    static Limits constant(double);  ///< Build dimension-0 limits (constant)

    bool operator<(const Limits&) const;         ///< Comparison operator
    Limits operator-() const;                    ///< Invert this limit
    Limits& operator=(const Limits&) = default;  ///< Assignment operator
    /// Equality operator
    inline bool operator==(const Limits& oth) const { return *this == static_cast<std::pair<double, double> >(oth); }
    inline bool operator!=(const Limits& oth) const { return !operator==(oth); }  ///< Inequality operator
    Limits& operator+=(double);                                                   ///< Add a constant to this limit
    Limits& operator-=(double);                                                   ///< Subtract a constant to this limit
    Limits& operator*=(double);                                                   ///< Multiply this limit by a constant
    friend Limits operator+(Limits, double);                                      ///< Add a constant to a limit
    friend Limits operator-(Limits, double);                                      ///< Subtract a constant to a limit
    friend Limits operator*(Limits, double);                                      ///< Multiply a limit by a constant

    Limits& validate();  ///< Ensure the limit object is valid by correcting it if necessary

    bool hasMin() const;                          ///< Have a lower limit?
    bool hasMax() const;                          ///< Have an upper limit?
    inline double min() const { return first; }   ///< Lower limit to apply on the variable
    inline double& min() { return first; }        ///< Lower limit to apply on the variable
    inline double max() const { return second; }  ///< Upper limit to apply on the variable
    inline double& max() { return second; }       ///< Upper limit to apply on the variable

    double x(double v) const;        ///< Find the [0,1] value scaled between minimum and maximum
    void in(double low, double up);  ///< Specify the lower and upper limits on the variable
    double range() const;            ///< Full variable range allowed

    Limits truncate(const Limits&) const;                              ///< Truncate limits to minimal/maximal values
    double trim(double) const;                                         ///< Limit a value to boundaries
    bool contains(double val, bool exclude_boundaries = false) const;  ///< Check if value is inside limits' boundaries
    Limits& apply(double (*)(double));                                 ///< Apply an operator on limits boundaries

    Limits compute(double (*)(double)) const;  ///< Compute a copy of limits with an operator applied on boundaries
    /// Compute a copy of limits with an operator applied on boundaries
    template <typename F>
    inline Limits compute(const F& op) const {
      return Limits{hasMin() ? op(min()) : INVALID, hasMax() ? op(max()) : INVALID};
    }
    bool valid() const;                                                    ///< Is there a lower and upper limit?
    inline const std::pair<double, double>& raw() const { return *this; }  ///< Raw value of limits

    /// Generate a collection of values from a number of bins
    /// \param[in] num_bins number of values to generate
    /// \param[in] log_scale generate according to a log10 scale?
    std::vector<double> generate(size_t num_bins, bool log_scale = false) const;
    /// Split the limits into sub-limits objects
    /// \param[in] num_bins number of sub-limits to generate
    /// \param[in] log_scale generate according to a log10 scale?
    std::vector<Limits> split(size_t num_bins, bool log_scale = false) const;

    friend std::ostream& operator<<(std::ostream&, const Limits&);  ///< Human-readable expression of the limits

    static constexpr double INVALID = -999.999;  ///< Invalid value placeholder (single-edged or invalid limits)
  };
}  // namespace cepgen

#endif
