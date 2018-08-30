#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/Event/Particle.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    incoming_beams( { { 6500., PDG::proton, 10 }, { 6500., PDG::proton, 10 } } ),
    mode( Mode::invalid ), structure_functions( new SF::SuriYennie )
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
    return incoming_beams.first.pz + incoming_beams.second.pz;
  }

  //-----------------------------------------------------------------------------------------------
  // User-friendly displayers
  //-----------------------------------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Mode& pm )
  {
    switch ( pm ) {
      case Kinematics::Mode::invalid:
        return os << "invalid";
      case Kinematics::Mode::ElectronElectron:
        return os << "electron/electron";
      case Kinematics::Mode::ElectronProton:
        return os << "electron/proton";
      case Kinematics::Mode::ProtonElectron:
        return os << "proton/electron";
      case Kinematics::Mode::ElasticElastic:
        return os << "elastic/elastic";
      case Kinematics::Mode::InelasticElastic:
        return os << "inelastic/elastic";
      case Kinematics::Mode::ElasticInelastic:
        return os << "elastic/inelastic";
      case Kinematics::Mode::InelasticInelastic:
        return os << "inelastic/inelastic";
    }
    return os;
  }

  //------------------------------------------------------------------------------------------------
  // List of kinematics limits
  //------------------------------------------------------------------------------------------------

  Kinematics::CutsList::CutsList()
  {
    initial.q2 = { 0., 1.e5 };
    central.pt_single.min() = 0.;
    remnants.mass_single = { 1.07, 320. };
  }
}

