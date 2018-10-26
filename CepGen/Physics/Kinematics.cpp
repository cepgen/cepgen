#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/KTFlux.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/Event/Particle.h"

#include <cmath>

namespace cepgen
{
  Kinematics::Kinematics() :
    incoming_beams( { { 6500., PDG::proton, KTFlux::invalid }, { 6500., PDG::proton, KTFlux::invalid } } ),
    mode( KinematicsMode::invalid ), structure_functions( new strfun::SuriYennie )
  {}

  Kinematics::~Kinematics()
  {}

  void
  Kinematics::setSqrtS( double sqrts )
  {
    incoming_beams.first.pz = incoming_beams.second.pz = 0.5 * sqrts;
  }

  double
  Kinematics::sqrtS() const
  {
    Particle::Momentum momA = Particle::Momentum::fromPxPyPzM( 0., 0., incoming_beams.first.pz, particleproperties::mass( incoming_beams.first.pdg ) );
    Particle::Momentum momB = Particle::Momentum::fromPxPyPzM( 0., 0., -incoming_beams.first.pz, particleproperties::mass( incoming_beams.second.pdg ) );
    std::cout << incoming_beams.second.pdg << "\t" << particleproperties::mass( incoming_beams.second.pdg ) << std::endl;
    const HeavyIon hi1 = (HeavyIon)incoming_beams.first.pdg, hi2 = (HeavyIon)incoming_beams.second.pdg;
    if ( hi1 || hi2 ) {
      momA *= hi1.A;
      momB *= hi2.A;
    }
      std::cout << momA << "\t" << momB << std::endl;
    return ( momA+momB ).mass();
  }

  //--------------------------------------------------------------------
  // User-friendly display of the kinematics mode
  //--------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const KinematicsMode& pm )
  {
    switch ( pm ) {
      case KinematicsMode::invalid:
        return os << "invalid";
      case KinematicsMode::ElectronElectron:
        return os << "electron/electron";
      case KinematicsMode::ElectronProton:
        return os << "electron/proton";
      case KinematicsMode::ProtonElectron:
        return os << "proton/electron";
      case KinematicsMode::ElasticElastic:
        return os << "elastic/elastic";
      case KinematicsMode::InelasticElastic:
        return os << "inelastic/elastic";
      case KinematicsMode::ElasticInelastic:
        return os << "elastic/inelastic";
      case KinematicsMode::InelasticInelastic:
        return os << "inelastic/inelastic";
    }
    return os;
  }

  //--------------------------------------------------------------------
  // User-friendly display of incoming particles
  //--------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Beam& beam )
  {
    if ( (HeavyIon)beam.pdg )
      os << (HeavyIon)beam.pdg;
    else
      os << beam.pdg;
    os << " (" << beam.pz << " GeV/c)";
    if ( beam.kt_flux != KTFlux::invalid )
      os << " [unint.flux: " << beam.kt_flux << "]";
    return os;
  }

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  Kinematics::CutsList::CutsList()
  {
    initial.q2 = { 0., 1.e5 };
    central.pt_single.min() = 0.;
    remnants.mass_single = { 1.07, 320. };
  }
}
