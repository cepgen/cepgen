#ifndef CepGen_Physics_Modes_h
#define CepGen_Physics_Modes_h

#include <iosfwd>

namespace cepgen
{
  /// Collection of enums for the definition of process mode
  namespace mode
  {
    /// Type of scattering
    enum class Kinematics
    {
      invalid            = 0,
      ElasticElastic     = 1, ///< proton-proton elastic case
      ElasticInelastic   = 2, ///< proton-proton single-dissociative (or inelastic-elastic) case
      InelasticElastic   = 3, ///< proton-proton single-dissociative (or elastic-inelastic) case
      InelasticInelastic = 4  ///< proton-proton double-dissociative case
    };

    /// Type of beam treatment
    enum class Beam
    {
      invalid          = 0,
      ProtonElastic    = 1, ///< Elastic scattering from proton
      ProtonInelastic  = 2, ///< Inelastic scattering from proton (according to the proton structure functions set)
      PointLikeScalar  = 3, ///< Trivial, spin-0 emission
      PointLikeFermion = 4, ///< Trivial, spin-1/2 emission
      CompositeScalar  = 5, ///< Composite pion emission
    };
  }
  /// Human-readable format of a process mode (elastic/dissociative parts)
  std::ostream& operator<<( std::ostream&, const mode::Kinematics& );
  /// Human-readable format of a beam mode (elastic/dissociative parts)
  std::ostream& operator<<( std::ostream&, const mode::Beam& );
}

#endif
