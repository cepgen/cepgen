#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

namespace cepgen
{
  enum class PDG;
  /// All useful properties about particles
  namespace particleproperties
  {
    /** \brief Mass of a particle, in GeV/c\f$^2\f$
     * \param pdg_id PDG identifier
     */
    double mass( const PDG& pdg_id );
    /** \brief Electric charge of a particle, in \f$e\f$
     * \param[in] pdg_id PDG id
     */
    double charge( const PDG& pdg_id );
    /** \brief Electric charge of a particle, in \f$e\f$
     * \param[in] id integer PDG id
     */
    double charge( int id );
    /** \brief Colour factor for a given particle
     * \param[in] pdg_id PDG id
     */
    unsigned short colours( const PDG& pdg_id );
    /** \brief Total decay width of an unstable particle, in GeV
     * \param[in] pdg_id PDG id
     */
    double width( const PDG& pdg_id );
  }
}

#endif
