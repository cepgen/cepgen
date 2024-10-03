/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2024  Laurent Forthomme
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
#include "CepGenPythia6/Pythia6Interface.h"

extern "C" {
double pyalem_(double& q2);
double pyalps_(double& q2);
extern double pymass_(int&);            ///< Get the particle's mass in GeV from Pythia
extern void pyexec_();                  ///< Launch the Pythia6 fragmentation
extern void pygive_(const char*, int);  ///< Set a parameter value to the Pythia module
extern void pyckbd_();
extern void pylist_(int&);              ///< List all the particles in the event in a human-readable format
extern void pyjoin_(int&, int&);        ///< Join two coloured particles in a colour singlet
extern void pyname_(int&, char*, int);  ///< Get a particle's human-readable name from Pythia
extern int pyk_(int&, int&);            ///< Get integer-valued event information from Pythia
extern double pyp_(int&, int&);         ///< Get real-valued event information from Pythia
extern int pychge_(int&);
void pystop_() { CG_INFO("pythia6:pystop") << "End of run"; }  ///< Purely virtual method to call at the end of the run
}

namespace cepgen::pythia6 {
  double pyalem(double q2) { return pyalem_(q2); }

  double pyalps(double q2) { return pyalps_(q2); }

  void pyexec() { pyexec_(); }

  int pychge(int pdgid) { return pychge_(pdgid); }

  void pyckbd() { pyckbd_(); }

  void pygive(const std::string& line) { pygive_(line.c_str(), line.length()); }

  void pyjoin(std::vector<int> join) {
    int njoin = join.size();
    return pyjoin_(njoin, *join.data());
  }

  int pyk(int id, int qty) { return pyk_(id, qty); }

  void pylist(int mlist) { pylist_(mlist); }

  double pymass(int pdgid) { return pymass_(pdgid); }

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

  double pyp(int id, int qty) { return pyp_(id, qty); }

  int pythia6Status(int cg_status) {
    switch (static_cast<cepgen::Particle::Status>(cg_status)) {
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

  int cepgenStatus(int py_status) {
    switch (py_status) {
      case 1:
        return static_cast<int>(cepgen::Particle::Status::FinalState);
      case 3:
        return static_cast<int>(cepgen::Particle::Status::Propagator);
      case 11:
        return static_cast<int>(cepgen::Particle::Status::Fragmented);
      case 21:
        return static_cast<int>(cepgen::Particle::Status::PrimordialIncoming);
      default:
        return py_status;
    }
  }

  void checkPDGid(int pdg_id) {
    if (cepgen::PDG::get().has(pdg_id))
      return;
    const auto name = pythia6::pyname(pdg_id);
    cepgen::ParticleProperties prop;
    prop.pdgid = pdg_id;
    prop.name = name;
    prop.human_name = name;
    //prop.colours = pyk(p + 1, 12);  // colour factor
    prop.mass = pymass(pdg_id);
    prop.width = -1.;  //pmas( pdg_id, 2 ),
    if (const auto ch = pychge(pdg_id); std::fabs(ch) > 0)
      prop.charges = {ch, -ch};
    prop.fermion = false;
    cepgen::PDG::get().define(prop);
  }
}  // namespace cepgen::pythia6
