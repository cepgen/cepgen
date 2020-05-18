#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Momentum.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Utils/String.h"

#include <cmath>

namespace cepgen
{
  Kinematics::Kinematics() :
    incoming_beams( { { 6500., PDG::proton, KTFlux::invalid }, { 6500., PDG::proton, KTFlux::invalid } } ),
    mode( KinematicsMode::invalid )
  {}

  Kinematics::Kinematics( const ParametersList& params ) :
    mode( (KinematicsMode)params.get<int>( "mode", (int)KinematicsMode::invalid ) )
  {
    incoming_beams.first.pdg = params.get<int>( "beam1id", 2212 );
    params.fill<double>( "beam1pz", incoming_beams.first.pz );
    const int hi_A1 = params.get<int>( "beam1A", 1 );
    const int hi_Z1 = params.get<int>( "beam1Z", 0 );
    if ( hi_Z1 != 0 )
      incoming_beams.first.pdg = HeavyIon( hi_A1, (Element)hi_Z1 );
    incoming_beams.second.pdg = params.get<int>( "beam2id", 2212 );
    params.fill<double>( "beam2pz", incoming_beams.second.pz );
    const int hi_A2 = params.get<int>( "beam2A", 1 );
    const int hi_Z2 = params.get<int>( "beam2Z", 0 );
    if ( hi_Z2 != 0 )
      incoming_beams.second.pdg = HeavyIon( hi_A2, (Element)hi_Z2 );
    const double sqrt_s = params.get<double>( "sqrts", -1. );
    if ( sqrt_s > 0. )
      setSqrtS( sqrt_s );
    for ( auto& lim : cuts.central.rawList() ) {
      params.fill<Limits>( lim.name, lim.limits );
      params.fill<double>( lim.name+"min", lim.limits.min() );
      params.fill<double>( lim.name+"max", lim.limits.max() );
    }
    //FIXME add the single-particles cuts parsing
    //--- outgoing remnants
    params.fill<Limits>( "mx", cuts.remnants.mass_single() );
    params.fill<double>( "mxmin", cuts.remnants.mass_single().min() );
    params.fill<double>( "mxmax", cuts.remnants.mass_single().max() );
    Limits xi_rng;
    params.fill<Limits>( "xi", xi_rng );
    if ( !xi_rng.valid() )
      xi_rng
        = { params.get<double>( "ximin", 0. ), params.get<double>( "ximax", 1. ) };
    if ( xi_rng.valid() )
      cuts.remnants.energy_single() = -( xi_rng-1. )*0.5*sqrtS();
  }

  void
  Kinematics::setSqrtS( double sqrts )
  {
    if ( incoming_beams.first.pdg != incoming_beams.second.pdg )
      throw CG_FATAL( "Kinematics" )
        << "Trying to set √s with asymmetric beams"
        << " (" << incoming_beams.first.pdg << "/" << incoming_beams.second.pdg << ").\n"
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
    return os;
  }

  //--------------------------------------------------------------------
  // List of kinematics limits
  //--------------------------------------------------------------------

  Kinematics::CutsList::CutsList()
  {
    initial.q2() = { 0., 1.e5 };
    central.pt_single().min() = 0.;
    remnants.mass_single() = { 1.07, 320. };
  }
}
