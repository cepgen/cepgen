/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2020-2025  Laurent Forthomme
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

// This file allows to link the MadGraph interfacing module without any process
// generation performed by MG5_aMC.
// Include it in your source file prior to any linking with libCepGenMadGraph.

#include "CepGen/Physics/Momentum.h"
#include "CepGenMadGraph/Process.h"
#include "CepGenMadGraph/ProcessFactory.h"

namespace cepgen::mg5amc {
  class DummyProcess final : public Process {
  public:
    using Process::Process;

    double eval() override { return 0.; }
    void initialise(const std::string&) override {}
    const std::vector<Momentum>& momenta() override { return momenta_; }

  private:
    std::vector<Momentum> momenta_;
  };
}  // namespace cepgen::mg5amc
REGISTER_MG5AMC_PROCESS("dummy", DummyProcess);
