/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024  Laurent Forthomme
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

#ifndef CepGen_Utils_ProcessVariablesAnalyser_h
#define CepGen_Utils_ProcessVariablesAnalyser_h

#include "CepGen/Core/SteeredObject.h"
#include "CepGen/Utils/Histogram.h"

namespace cepgen {
  namespace proc {
    class Process;
  }
  namespace utils {
    class ProcessVariablesAnalyser : public SteeredObject<ProcessVariablesAnalyser> {
    public:
      static ProcessVariablesAnalyser& get(const ParametersList& = ParametersList());
      static ParametersDescription description();

      void reset(const proc::Process&);
      void analyse(proc::Process&, double weight);

    private:
      explicit ProcessVariablesAnalyser(const ParametersList&);
      std::unordered_map<std::string, Hist1D> hists_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
