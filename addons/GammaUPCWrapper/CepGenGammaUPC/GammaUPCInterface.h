/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2024-2026  Laurent Forthomme
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

#ifndef CepGenGammaUPC_GammaUPCInterface_h
#define CepGenGammaUPC_GammaUPCInterface_h

#include <array>
#include <iosfwd>

namespace gammaUPC {
  struct to_collider_t {
    std::array<double, 2> ebeam, xbk, q2fact;
    std::array<int, 2> lpp;
  };
  struct to_heavyion_pdg_t {
    std::array<int, 2> nb_proton, nb_neutron;
  };
}  // namespace gammaUPC

extern "C" {
extern gammaUPC::to_collider_t to_collider_;
}

namespace gammaUPC {
  enum struct BeamsContent { none = 0, heavyion_cb, module_var };
  std::ostream& operator<<(std::ostream&, const BeamsContent&);
  BeamsContent setBeamsContent(const std::array<int, 2>& z, const std::array<int, 2>& a);
}  // namespace gammaUPC

namespace gammaUPC {
  enum struct HeavyIonMode { HardSphere, WoodsSaxon };
  double photonFluxP(double x, double gamma);
  double photonFluxA(double a, double gamma, double z, double ra);
  double twoPhotonFluxPP(double x1, double x2, int force_p_no_hadronisation);
  double twoPhotonFluxPA(double x1, double x2, int force_p_no_hadronisation, HeavyIonMode mode);
  double twoPhotonFluxAB(double x1, double x2, int force_p_no_hadronisation, HeavyIonMode mode);
}  // namespace gammaUPC

#endif
