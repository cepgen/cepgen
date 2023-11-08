/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2023  Laurent Forthomme
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

#include "CepGen/Core/Exception.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/Pythia6Wrapper/Pythia6Interface.h"

extern "C" {
/// Get the particle's mass in GeV from Pythia
extern double pymass_(int&);
/// Launch the Pythia6 fragmentation
extern void pyexec_();
/// Set a parameter value to the Pythia module
extern void pygive_(const char*, int);
extern void pyckbd_();
/// List all the particles in the event in a human-readable format
extern void pylist_(int&);
/// Join two coloured particles in a colour singlet
extern void pyjoin_(int&, int&);
/// Get a particle's human-readable name from Pythia
extern void pyname_(int&, char*, int);
/// Get integer-valued event information from Pythia
extern int pyk_(int&, int&);
/// Get real-valued event information from Pythia
extern double pyp_(int&, int&);
extern double pychge_(int&);
/// Purely virtual method to call at the end of the run
void pystop_() { CG_INFO("pythia6:pystop") << "End of run"; }
}

namespace pythia6 {
  void pyexec() { pyexec_(); }

  double pymass(int pdgid) { return pymass_(pdgid); }

  void pyckbd() { pyckbd_(); }

  void pygive(const std::string& line) { pygive_(line.c_str(), line.length()); }

  void pylist(int mlist) { pylist_(mlist); }

  int pyk(int id, int qty) { return pyk_(id, qty); }

  double pyp(int id, int qty) { return pyp_(id, qty); }

  double pychge(int pdgid) { return pychge_(pdgid); }

  std::string pyname(int pdgid) {
    // maximal number of characters to fetch for the particle's name
    static constexpr unsigned short NAME_CHR = 16;

    char out[NAME_CHR];
    std::string s;
    pyname_(pdgid, out, NAME_CHR);
    s = std::string(out, NAME_CHR);
    s.erase(remove(s.begin(), s.end(), ' '), s.end());
    return s;
  }

  void pyjoin(std::vector<int> join) {
    int njoin = join.size();
    return pyjoin_(njoin, *join.data());
  }

  int pythia6Status(int cg_status) {
    switch ((cepgen::Particle::Status)cg_status) {
      case cepgen::Particle::Status::PrimordialIncoming:
        return 21;
      case cepgen::Particle::Status::FinalState:
      case cepgen::Particle::Status::Undecayed:
        return 1;
      case cepgen::Particle::Status::Unfragmented:
        return 3;
      case cepgen::Particle::Status::Fragmented:
      case cepgen::Particle::Status::Propagator:
      case cepgen::Particle::Status::Incoming:
        return 11;
      default:
        throw CG_FATAL("pythia6:status") << "No conversion rule for CepGen status code: " << cg_status << ".";
    }
  }

  cepgen::Particle::Status cepgenStatus(int py_status) {
    CG_LOG << py_status;
    switch (py_status) {
      case 1:
        return cepgen::Particle::Status::FinalState;
      case 3:
        return cepgen::Particle::Status::Propagator;
      case 11:
        return cepgen::Particle::Status::Fragmented;
      case 21:
        return cepgen::Particle::Status::PrimordialIncoming;
      default:
        return (cepgen::Particle::Status)py_status;
    }
  }

  void checkPDGid(int pdg_id) {
    if (cepgen::PDG::get().has(pdg_id))
      return;
    const auto name = pythia6::pyname(pdg_id);
    cepgen::ParticleProperties prop;
    prop.pdgid = pdg_id;
    prop.name = name;
    prop.descr = name;
    //prop.colours = pyk(p + 1, 12);  // colour factor
    prop.mass = pymass(pdg_id);
    prop.width = -1.;              //pmas( pdg_id, 2 ),
    prop.charge = pychge(pdg_id);  // charge
    prop.fermion = false;
    cepgen::PDG::get().define(prop);
  }
}  // namespace pythia6
