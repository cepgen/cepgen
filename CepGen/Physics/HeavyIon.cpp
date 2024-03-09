/*
 *  CepGen: a central exclusive processes event generator
 *  Copyright (C) 2013-2022  Laurent Forthomme
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

#include <math.h>

#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

#define ELEM_STR(x) #x
#define ELEM_EVAL(x) ELEM_STR(x)
#define DEF_ELEM(x) \
  case Element::x:  \
    return os << ELEM_EVAL(x)

namespace cepgen {
  HeavyIon::HeavyIon(unsigned short a, const Element& z) : Z(z), A(a) {}

  HeavyIon HeavyIon::fromPdgId(pdgid_t pdg) {
    if (pdg == PDG::neutron)
      return neutron();
    else if (pdg == PDG::proton)
      return proton();
    else if (pdg / 10'000'000 != 0)
      return HeavyIon(pdg % 1000, static_cast<Element>((pdg / 1000) % 1000));
    CG_WARNING("HeavyIon") << "Failed to parse heavy ion from PDG id=" << pdg << ".";
    return HeavyIon(0, Element::invalid);
  }

  HeavyIon::operator pdgid_t() const {
    // Pythia8 convention/10-1e10+1e6
    if (*this == proton())
      return PDG::proton;
    if (*this == neutron())
      return PDG::neutron;
    return (pdgid_t)(10'000'000 + 1000 * (unsigned short)Z + A);
  }

  bool HeavyIon::isHI(const spdgid_t& pdgid) { return pdgid / 10'000'000 != 0; }

  bool HeavyIon::isHI(const ParticleProperties& prop) {
    return isHI(prop.pdgid);  //FIXME can refine a bit
  }

  double HeavyIon::massP() const {
    if (Z == Element::invalid)
      throw CG_FATAL("HeavyIon:massP") << "Invalid heavy ion: " << (*this) << "!";
    return (short)Z * PDG::get().mass(PDG::proton);
  }

  double HeavyIon::massN() const {
    if (Z == Element::invalid)
      throw CG_FATAL("HeavyIon:massN") << "Invalid heavy ion: " << (*this) << "!";
    return (A - (short)Z) * PDG::get().mass(PDG::neutron);
  }

  double HeavyIon::radius() const {
    switch (Z) {
      case Element::H: {
        if (A == 1)  // simple proton
          return 0.841e-15;
        else  // deuteron
          return 2.128e-15;
      }
      case Element::Cu:
        return 4.214e-15;
      case Element::Xe:
        return 5.36e-15;
      case Element::Au:
        return 6.38e-15;
      case Element::Pb:
        return 6.624e-15;
      default: {
        CG_WARNING("HeavyIon:radius") << "Using hard-sphere approximation R ~ 1.2 A^(1/3).";
        return 1.2e-15 * cbrt(A);
      }
    }
  }

  std::ostream& operator<<(std::ostream& os, const HeavyIon& hi) {
    if (hi == HeavyIon::proton())
      return os << "proton";
    if (hi == HeavyIon::neutron())
      return os << "neutron";
    std::ostringstream oss;
    oss << hi.Z;
    if (oss.str().empty() || hi.Z == Element::invalid)
      return os << "HI{Z=" << (unsigned short)hi.Z << ", A=" << hi.A << "}";
    return os << hi.A << oss.str();
  }

  std::ostream& operator<<(std::ostream& os, const Element& elem) {
    switch (elem) {
      DEF_ELEM(invalid);
      DEF_ELEM(neutron);
      DEF_ELEM(H);
      DEF_ELEM(C);
      DEF_ELEM(O);
      DEF_ELEM(Al);
      DEF_ELEM(Cu);
      DEF_ELEM(Xe);
      DEF_ELEM(Au);
      DEF_ELEM(Pb);
      DEF_ELEM(U);
    }
    return os;
  }
}  // namespace cepgen

#undef DEF_ELEM
#undef ELEM_STR
#undef ELEM_EVAL
