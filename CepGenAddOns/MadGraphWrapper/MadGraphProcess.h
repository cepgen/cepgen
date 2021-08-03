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

#ifndef CepGenAddOns_MadGraphWrapper_MadGraphProcess_h
#define CepGenAddOns_MadGraphWrapper_MadGraphProcess_h

#include <memory>
#include <string>

#include "CepGen/Physics/Momentum.h"

// forward-declaration of base MadGraph standalone_cpp process
class CPPProcess;

namespace cepgen {
  /// Wrapper around a generic MadGraph CPPProcess definition
  class MadGraphProcess {
  public:
    MadGraphProcess();
    ~MadGraphProcess();

    const std::string& name() const { return name_; }
    const std::string& description() const { return descr_; }

    void initialise(const std::string&);
    double eval();

    MadGraphProcess& setMomentum(size_t, const Momentum&);
    const std::vector<Momentum>& momenta();
    const std::vector<double>& masses() const;

  private:
    std::unique_ptr<CPPProcess> proc_;
    std::vector<double*> mom_;

    const std::string name_;
    const std::string descr_;
    const std::array<int, 2> incoming_pdgids_;
    const std::vector<int> central_pdgids_;

    std::vector<Momentum> momenta_;
  };
}  // namespace cepgen

#endif
