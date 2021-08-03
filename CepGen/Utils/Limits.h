/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2021  Laurent Forthomme
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

namespace cepgen {
  /// Validity interval for a variable
  class Limits : private std::pair<double, double> {
  public:
    /// Define lower and upper limits on a quantity
    Limits(double min = INVALID, double max = INVALID);
    /// Copy constructor
    Limits(const Limits&);

    Limits operator-() const;                       ///< Invert this limit
    Limits& operator=(const Limits&) = default;     ///< Assignment operator
    Limits& operator+=(double c);                   ///< Add a constant to this limit
    Limits& operator-=(double c);                   ///< Subtract a constant to this limit
    Limits& operator*=(double c);                   ///< Multiply this limit by a constant
    friend Limits operator+(Limits lim, double c);  ///< Add a constant to a limit
    friend Limits operator-(Limits lim, double c);  ///< Subtract a constant to a limit
    friend Limits operator*(Limits lim, double c);  ///< Multiply a limit by a constant

    /// Ensure the limit object is valid by correcting it if necessary
    Limits& validate();

    /// Lower limit to apply on the variable
    double min() const { return first; }
    /// Lower limit to apply on the variable
    double& min() { return first; }
    /// Upper limit to apply on the variable
    double max() const { return second; }
    /// Upper limit to apply on the variable
    double& max() { return second; }
    /// Export the limits into external variables
    void save(bool& on, double& lmin, double& lmax) const;
    /// Find the [0,1] value scaled between minimum and maximum
    double x(double v) const;
    /// Specify the lower and upper limits on the variable
    void in(double low, double up);
    /// Full variable range allowed
    double range() const;
    /// Have a lower limit?
    bool hasMin() const;
    /// Have an upper limit?
    bool hasMax() const;
    /// Check if the value is inside limits' boundaries
    bool contains(double val) const;
    /// Is there a lower and upper limit?
    bool valid() const;
    /// Raw value of the limits
    const std::pair<double, double> raw() const { return *this; }

    /// Human-readable expression of the limits
    friend std::ostream& operator<<(std::ostream&, const Limits&);

    /// Placeholder for an invalid value in a limit (for single-edged or invalid limits)
    static constexpr double INVALID = -999.999;
  };
}  // namespace cepgen

#endif
