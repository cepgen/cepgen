/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2025  Laurent Forthomme
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

#include "CepGen/Process/FactorisedProcess.h"
#include "CepGenMadGraph/Process.h"

namespace cepgen::mg5amc {
  class ProcessBuilder : public proc::FactorisedProcess {
  public:
    explicit ProcessBuilder(const ParametersList& params, bool load_library = true);

    static ParametersDescription description();

    void addEventContent() override;

  protected:
    void loadMG5Library() const;
    void prepareSteeringCard() const;

    std::unique_ptr<mg5amc::Process> mg5_proc_;
  };
}  // namespace cepgen::mg5amc
