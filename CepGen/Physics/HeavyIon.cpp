#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"

#define DEF_ELEM(x) \
  case Element::x:  \
    return os << "x"

namespace cepgen {
  HeavyIon::HeavyIon(pdgid_t pdg)
      : Z(pdg / 1000000 == 0 ? Element::invalid : (Element)((pdg / 1000) % 1000)),
        A((Z != Element::invalid) ? pdg % 1000 : 0) {}

  HeavyIon::operator pdgid_t() const {
    // Pythia8 convention/10-1e10+1e6
    return (pdgid_t)(1000000 + 1000 * (unsigned short)Z + A);
  }

  HeavyIon::operator bool() const {
    return Z != Element::invalid;  // skip the proton
  }

  double HeavyIon::mass(const HeavyIon& hi) {
    if (!hi)
      throw CG_FATAL("mass") << "Invalid heavy ion: " << hi << "!";
    return (short)hi.Z * PDG::get().mass(PDG::proton);
  }

  std::ostream& operator<<(std::ostream& os, const HeavyIon& hi) {
    std::ostringstream oss;
    oss << hi.Z;
    if (oss.str().empty())
      return os << "HI{Z=" << (unsigned short)hi.Z << ", A=" << hi.A << "}";
    return os << hi.A << oss.str();
  }

  std::ostream& operator<<(std::ostream& os, const Element& elem) {
    switch (elem) {
      DEF_ELEM(invalid);
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
