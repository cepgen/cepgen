#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"

#include "CepGen/StructureFunctions/SuriYennie.h"
#include "CepGen/Event/Particle.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    inp( { 6500., 6500. } ), inpdg( { PDG::Proton, PDG::Proton } ),
    mode( Mode::ElasticElastic ), structure_functions( new SF::SuriYennie )
  {}

  /*Kinematics::Kinematics( const Kinematics& kin ) :
    inp( kin.inp ), inpdg( kin.inpdg ), mode( kin.mode ),
    structure_functions( kin.structure_functions.get() )
  {}*/

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

