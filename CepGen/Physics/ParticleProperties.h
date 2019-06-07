#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

namespace cepgen
{
  enum class PDG;
  struct HeavyIon;
  /// All useful properties about particles
  namespace particleproperties
  {
    /// Mass of a heavy ion, in GeV/c\f$^2\f$
    /// \param hi Heavy ion type
    double mass( const HeavyIon& hi );
    bool isFermion( const PDG& pdg_id );
  }
}

#endif
