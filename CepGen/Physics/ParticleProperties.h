#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

#include <string>
#include <iosfwd>

namespace cepgen {
  /// Alias for the integer-like particle PDG id
  typedef unsigned long pdgid_t;
  /// A collection of physics constants associated to a single particle
  struct ParticleProperties {
    pdgid_t pdgid;            ///< PDG identifier
    std::string name;         ///< Particle name
    std::string description;  ///< Human-readable name
    short colours;            ///< Colour factor
    double mass;              ///< Mass, in GeV/c\f$^2\f$
    double width;             ///< Decay width, in GeV/c\f$^2\f$
    short charge;             ///< Electric charge, in \f$e\f$/3
    bool fermion;             ///< Is the particle a fermion?
    friend std::ostream& operator<<(std::ostream&, const ParticleProperties&);
  };
}  // namespace cepgen

#endif
