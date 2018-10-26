#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

namespace cepgen
{
  enum class PDG;
  struct HeavyIon;
  /// All useful properties about particles
  namespace particleproperties
  {
    /// Mass of a particle, in GeV/c\f$^2\f$
    /// \param pdg_id PDG identifier
    double mass( const PDG& pdg_id );
    /// Mass of a heavy ion, in GeV/c\f$^2\f$
    /// \param hi Heavy ion type
    double mass( const HeavyIon& hi );
    /// Electric charge of a particle, in \f$e\f$
    /// \param[in] pdg_id PDG id
    double charge( const PDG& pdg_id );
    /// Electric charge of a particle, in \f$e\f$
    /// \param[in] id integer PDG id
    double charge( int id );
    /// Colour factor for a given particle
    /// \param[in] pdg_id PDG id
    unsigned short colours( const PDG& pdg_id );
    /// Total decay width of an unstable particle, in GeV
    /// \param[in] pdg_id PDG id
    double width( const PDG& pdg_id );
  }
}

#endif
