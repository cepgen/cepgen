#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

namespace CepGen
{
  enum class PDG;
  namespace ParticleProperties
  {
    /// Mass (in GeV) of a particle
    /// \param pdg_id PDG identifier
    /// \return Mass of the particle in GeV/c\f${}^2\f$
    double mass( const PDG& pdg_id );
    /// Electric charge of a particle, in \f$e\f$
    /// \param[in] pdg_id PDG id
    double charge( const PDG& pdg_id );
    /// Electric charge of a particle, in \f$e\f$
    /// \param[in] id integer PDG id
    double charge( int id );
    unsigned short colours( const PDG& pdg_id );
    /// Total decay width of an unstable particle, in GeV
    /// \param[in] pdg_id PDG (PDG ID)
    double width( const PDG& pdg_id );
  }
}

#endif

