#include "Kinematics.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    inp( { 6500., 6500. } ), inpdg( { PDG::Proton, PDG::Proton } ),
    inhi( { HeavyIon::Proton(), HeavyIon::Proton() } ), kt_fluxes( { 10, 10 } ),
    mode( Mode::ElasticElastic ), structure_functions( StructureFunctions::SuriYennie )
  {}

  Kinematics::~Kinematics()
  {}

  void
  Kinematics::setSqrtS( double sqrts )
  {
    const double pin = 0.5 * sqrts;
    inp = { pin, pin };
  }

  double
  Kinematics::sqrtS() const
  {
    return inp.first + inp.second;
  }

  //-----------------------------------------------------------------------------------------------
  // User-friendly displayers
  //-----------------------------------------------------------------------------------------------

  std::ostream&
  operator<<( std::ostream& os, const Kinematics::Mode& pm )
  {
    switch ( pm ) {
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

