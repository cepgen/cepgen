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

#ifndef CepGenAddOns_Pythia6Wrapper_Pythia6Interface_h
#define CepGenAddOns_Pythia6Wrapper_Pythia6Interface_h

#include <string>
#include <vector>

extern "C" {
/// Particles content of the event
struct pyjets_t {
  int n;  ///< Number of particles in the event
  int npad;
  int k[5][4000];     ///< Particles' general information (status, PDG id, mother, daughter 1, daughter 2)
  double p[5][4000];  ///< Particles' kinematics, in GeV (px, py, pz, E, M)
  double v[5][4000];  ///< Primary vertex location for the particles
};
extern pyjets_t pyjets_;
}

/// Pythia 6 utilities namespace
namespace pythia6 {
  int pythia6Status(int);
  int cepgenStatus(int);
  void checkPDGid(int);

  void pyexec();
  double pymass(int pdgid_);
  void pyckbd();
  void pygive(const std::string&);
  void pylist(int mlist);
  int pyk(int id, int qty);
  double pyp(int id, int qty);
  std::string pyname(int pdgid);
  /// Connect entries with colour flow information
  /// \param[in] join List of particles unique identifier to join in the colour flow
  void pyjoin(std::vector<int> join);
}  // namespace pythia6

#endif
