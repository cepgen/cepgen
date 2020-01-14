#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/Momentum.h"

#include "CepGen/StructureFunctions/Parameterisation.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include <cmath>

namespace cepgen
{
  Kinematics::Kinematics() :
    mode( KinematicsMode::invalid )
  {}

  void
  Kinematics::setSqrtS( double sqrts )
  {
    if ( incoming_beams.first.pdg != incoming_beams.second.pdg )
      throw CG_FATAL( "Kinematics" )
        << "Trying to set √s with asymmetric beams.\n"
        << "Please fill incoming beams objects manually!";
    incoming_beams.first.pz = incoming_beams.second.pz = 0.5 * sqrts;
  }

  double
  Kinematics::sqrtS() const
  {
    const HeavyIon hi1( incoming_beams.first.pdg ), hi2( incoming_beams.second.pdg );
    const double m1 = hi1 ? HeavyIon::mass( hi1 ) : PDG::get().mass( incoming_beams.first .pdg );
    const double m2 = hi2 ? HeavyIon::mass( hi2 ) : PDG::get().mass( incoming_beams.second.pdg );
    const auto p1 = Momentum::fromPxPyPzM( 0., 0., +incoming_beams.first .pz, m1 );
    const auto p2 = Momentum::fromPxPyPzM( 0., 0., -incoming_beams.second.pz, m2 );
    return ( p1+p2 ).mass();
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
      os << PDG::get().name( beam.pdg );
    os << " (" << beam.pz << " GeV/c)";
    if ( beam.kt_flux != KTFlux::invalid )
      os << " [unint.flux: " << beam.kt_flux << "]";
    if ( beam.form_factors )
      os << ", " << *beam.form_factors;
    return os;
  }

  //--------------------------------------------------------------------
  // Default beam parameters
  //--------------------------------------------------------------------

  Kinematics::Beam::Beam() :
    pz( 0. ), pdg( PDG::proton ), kt_flux( KTFlux::invalid ),
    form_factors( new ff::StandardDipole )
  {}

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
