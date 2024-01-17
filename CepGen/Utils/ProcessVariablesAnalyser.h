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
    class Drawer;
    class ProcessVariablesAnalyser : public SteeredObject<ProcessVariablesAnalyser> {
    public:
      explicit ProcessVariablesAnalyser(const proc::Process&, const ParametersList&);

      static ParametersDescription description();

      void feed(double weight);
      void analyse();

    private:
      const proc::Process& proc_;
      const std::unique_ptr<Drawer> drawer_;
      std::unordered_map<std::string, Hist1D> hists_;
    };
  }  // namespace utils
}  // namespace cepgen

#endif
