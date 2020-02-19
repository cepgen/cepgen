#ifndef CepGen_Physics_KinematicsMode_h
#define CepGen_Physics_KinematicsMode_h

#include <iosfwd>

namespace cepgen
{
  /// Type of kinematics to consider for the process
  enum class KinematicsMode
  {
    invalid = -1,
    ElectronProton = 0,     ///< electron-proton elastic case
    ElasticElastic = 1,     ///< proton-proton elastic case
    ElasticInelastic = 2,   ///< proton-proton single-dissociative (or inelastic-elastic) case
    InelasticElastic = 3,   ///< proton-proton single-dissociative (or elastic-inelastic) case
    InelasticInelastic = 4, ///< proton-proton double-dissociative case
    ProtonElectron,
    ElectronElectron
  };
  /// Human-readable format of a process mode (elastic/dissociative parts)
  std::ostream& operator<<( std::ostream&, const KinematicsMode& );
}

#endif

