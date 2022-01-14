/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#ifndef CepGen_Utils_Histogram_h
#define CepGen_Utils_Histogram_h

#include <string>

namespace cepgen {
  namespace utils {
    /**
     * \brief Generic text-based plotting utility
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class Histogram {
    public:
      Histogram() = default;
      virtual ~Histogram() = default;

      /// Reset the histogram
      virtual void clear() = 0;
      /// Rescale all histogram bins by a constant factor
      virtual void scale(double) = 0;
      /// Compute the histogram integral
      virtual double integral() const = 0;
      /// Retrieve the maximum bin value
      virtual double minimum() const = 0;
      /// Retrieve the minimum bin value
      virtual double maximum() const = 0;
      /// Normalise the histogram to a given constant
      void normalise(double integ = 1.) { scale(integ / integral()); }
    };
  }  // namespace utils
}  // namespace cepgen

#endif
