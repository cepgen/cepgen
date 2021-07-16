#include <iostream>

#include "CepGen/Physics/Modes.h"

namespace cepgen {
  namespace mode {
    std::ostream& operator<<(std::ostream& os, const Kinematics& pm) {
      switch (pm) {
        case Kinematics::invalid:
          return os << "{invalid}";
        case Kinematics::ElasticElastic:
          return os << "elastic/elastic";
        case Kinematics::InelasticElastic:
          return os << "inelastic/elastic";
        case Kinematics::ElasticInelastic:
          return os << "elastic/inelastic";
        case Kinematics::InelasticInelastic:
          return os << "inelastic/inelastic";
      }
      return os;
    }

    std::ostream& operator<<(std::ostream& os, const Beam& type) {
      switch (type) {
        case Beam::invalid:
          return os << "{invalid}";
        case Beam::ProtonElastic:
          return os << "el.proton";
        case Beam::PointLikeScalar:
          return os << "gen.scalar";
        case Beam::PointLikeFermion:
          return os << "gen.fermion";
        case Beam::CompositeScalar:
          return os << "comp.scalar";
        case Beam::ProtonInelastic:
          return os << "inel.proton";
      }
      return os;
    }
  }  // namespace mode
}  // namespace cepgen
