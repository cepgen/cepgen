#include "Kinematics.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/utils.h"
#include "CepGen/Event/Particle.h"
#include "CepGen/Physics/PDG.h"

namespace CepGen
{
  Kinematics::Kinematics() :
    inp( { 6500., 6500. } ), inpdg( { PDG::Proton, PDG::Proton } ),
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

  void
  Kinematics::dump( std::ostream& os ) const
  {
    os << std::setfill(' ');
    os << "===== Central system\n";
    for ( const auto& pdg_lim : cuts.central )
      os << std::setw(30) << pdg_lim.first << ": " << pdg_lim.second;

    os << "===== Initial state\n";
    for ( const auto& pdg_lim : cuts.initial )
      os << std::setw(30) << pdg_lim.first << ": " << pdg_lim.second;

    os << "===== Remnants\n";
    for ( const auto& pdg_lim : cuts.remnants )
      os << std::setw(30) << pdg_lim.first << ": " << pdg_lim.second;
  }

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

  Kinematics::CutsList::CutsList() :
    initial( { { Cuts::q2, { 0., 1.e5 } } } ),
    central( { { Cuts::pt_single, 0. } } ),
    remnants( { { Cuts::mass_single, { 1.07, 320. } } } )
  {}
}

