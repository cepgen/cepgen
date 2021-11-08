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

#ifndef CepGen_Integration_Integrand_h
#define CepGen_Integration_Integrand_h

#include <functional>
#include <memory>
#include <vector>

namespace cepgen {
  class Parameters;
  namespace utils {
    class Timer;
  }
  /// Wrapper to the function to be integrated
  class Integrand {
  public:
    explicit Integrand(const Parameters*);
    virtual ~Integrand();

    void setFunction(size_t ndim, const std::function<double(const std::vector<double>&)>& func);

    /// Compute the integrand for a given phase space point
    virtual double eval(const std::vector<double>& x);
    /// Phase space dimension
    virtual size_t size() const { return gen_integr_.ndim; }

  protected:
    const Parameters* params_;                 ///< Generator-owned runtime parameters
    const std::unique_ptr<utils::Timer> tmr_;  ///< A precious timekeeper for event timing

  private:
    struct GenericIntegrand {
      std::function<double(const std::vector<double>&)> function;
      size_t ndim;
    } gen_integr_;
  };
}  // namespace cepgen

#endif
