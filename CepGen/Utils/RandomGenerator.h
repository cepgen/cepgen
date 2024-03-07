/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2023  Laurent Forthomme
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

#ifndef CepGen_Utils_RandomGenerator_h
#define CepGen_Utils_RandomGenerator_h

#include <array>
#include <string>
#include <vector>

#include "CepGen/Core/SteeredObject.h"

namespace cepgen {
  namespace utils {
    /// A random number generator
    /// \author L. Forthomme <laurent.forthomme@cern.ch>
    /// \date Nov 2023
    class RandomGenerator : public SteeredObject<RandomGenerator> {
    public:
      /// Default constructor
      explicit RandomGenerator(const ParametersList&);

      static ParametersDescription description();

      // base distributions
      virtual int uniformInt(int min, int max) = 0;
      virtual double uniform(double min = 0., double max = 1.) = 0;
      virtual double normal(double mean = 0., double rms = 1.) = 0;

      // specialised distributions
      virtual double exponential(double exponent = 1.);
      virtual double breitWigner(double mean = 0., double scale = 1.);
      virtual double landau(double location = 0., double width = 1.);
      virtual int poisson(double mean = 0.);

      /// Retrieve the engine object
      template <typename T>
      T* engine() {
        return static_cast<T*>(enginePtr());
      }

    protected:
      unsigned long long seed_;
      virtual void* enginePtr();  ///< engine object
    };
  }  // namespace utils
}  // namespace cepgen

#endif
