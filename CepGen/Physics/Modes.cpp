#include "CepGen/Physics/Modes.h"

#include <iostream>

namespace cepgen {
  std::ostream& operator<<(std::ostream& os, const mode::Kinematics& pm) {
    switch (pm) {
      case mode::Kinematics::invalid:
        return os << "{invalid}";
      case mode::Kinematics::ElasticElastic:
        return os << "elastic/elastic";
      case mode::Kinematics::InelasticElastic:
        return os << "inelastic/elastic";
      case mode::Kinematics::ElasticInelastic:
        return os << "elastic/inelastic";
      case mode::Kinematics::InelasticInelastic:
        return os << "inelastic/inelastic";
    }
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const mode::Beam& type) {
    switch (type) {
      case mode::Beam::invalid:
        return os << "{invalid}";
      case mode::Beam::Electron:
        return os << "electron";
      case mode::Beam::ProtonElastic:
        return os << "el.proton";
      case mode::Beam::PointLikeScalar:
        return os << "gen.scalar";
      case mode::Beam::PointLikeFermion:
        return os << "gen.fermion";
      case mode::Beam::CompositeScalar:
        return os << "comp.scalar";
      case mode::Beam::ProtonInelastic:
        return os << "inel.proton";
    }
    return os;
  }
}  // namespace cepgen
