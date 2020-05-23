#include "CepGen/Physics/Kinematics.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/HeavyIon.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/Momentum.h"

#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"

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
    //--- incoming beams
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
    //--- structure functions
    auto strfun = params.get<ParametersList>( "structureFunctions" );
    if ( !strfun.empty() || !structure_functions ) {
      if ( strfun.name<int>( -999 ) == -999 )
        strfun.setName<int>( 11 ); // default is Suri-Yennie
      structure_functions = strfun::StructureFunctionsFactory::get().build( strfun );
    }
    //--- phase space definition
    for ( auto& lim : cuts.central.rawList() ) {
      params.fill<Limits>( lim.name, lim.limits );
      params.fill<double>( lim.name+"min", lim.limits.min() );
      params.fill<double>( lim.name+"max", lim.limits.max() );
    }
    //--- parton fluxes for kt-factorisation
    if ( params.has<std::vector<int> >( "ktFluxes" ) ) {
      auto kt_fluxes = params.get<std::vector<int> >( "ktFluxes" );
      if ( !kt_fluxes.empty() ) {
        incoming_beams.first.kt_flux = (KTFlux)kt_fluxes.at( 0 );
        incoming_beams.second.kt_flux = ( kt_fluxes.size() > 1 )
          ? (KTFlux)kt_fluxes.at( 1 )
          : (KTFlux)kt_fluxes.at( 0 );
      }
    }
    else if ( params.has<int>( "ktFluxes" ) )
      incoming_beams.first.kt_flux = incoming_beams.second.kt_flux
        = (KTFlux)params.get<int>( "ktFluxes" );
    //FIXME add the single-particles cuts parsing
    //--- outgoing remnants
    params.fill<Limits>( "mx", cuts.remnants.mass_single() );
    params.fill<double>( "mxmin", cuts.remnants.mass_single().min() );
    params.fill<double>( "mxmax", cuts.remnants.mass_single().max() );
    params.fill<Limits>( "yj", cuts.remnants.rapidity_single() );
    params.fill<double>( "yjmin", cuts.remnants.rapidity_single().min() );
    params.fill<double>( "yjmax", cuts.remnants.rapidity_single().max() );
    Limits xi_rng;
    params.fill<Limits>( "xi", xi_rng );
    if ( !xi_rng.valid() )
      xi_rng
        = { params.get<double>( "ximin", 0. ), params.get<double>( "ximax", 1. ) };
    if ( xi_rng.valid() ) {
      if ( incoming_beams.first.pz != incoming_beams.second.pz )
        CG_WARNING( "Kinematics" )
          << "xi range setting is not yet supported for asymmetric beams.\n\t"
          << "Handle with care.";
      cuts.remnants.energy_single() = -( xi_rng-1. )*0.5*sqrtS();
    }
  }

  Kinematics&
  Kinematics::setSqrtS( double sqrts )
  {
    if ( incoming_beams.first.pdg != incoming_beams.second.pdg )
      throw CG_FATAL( "Kinematics" )
        << "Trying to set √s with asymmetric beams"
        << " (" << incoming_beams.first.pdg << "/" << incoming_beams.second.pdg << ").\n"
        << "Please fill incoming beams objects manually!";
    incoming_beams.first.pz = incoming_beams.second.pz = 0.5 * sqrts;
    return *this;
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

  Kinematics&
  Kinematics::setStructureFunctions( int sf_model, int sr_model )
  {
    auto sf_params = ParametersList()
      .setName<int>( sf_model )
      .set<ParametersList>( "sigmaRatio", ParametersList()
        .setName<int>( sr_model ) );
    const unsigned long kLHAPDFCodeDec = 10000000, kLHAPDFPartDec = 1000000;
    if ( sf_model / kLHAPDFCodeDec == 1 ) { // SF from parton
      const unsigned long icode = sf_model % kLHAPDFCodeDec;
      sf_params
        .setName<int>( (int)strfun::Type::Partonic )
        .set<int>( "pdfId", icode % kLHAPDFPartDec )
        .set<int>( "mode", icode / kLHAPDFPartDec ); // 0, 1, 2
    }
    structure_functions = strfun::StructureFunctionsFactory::get().build( sf_params );
    return *this;
  }

  ParametersList
  Kinematics::parameters() const
  {
    ParametersList params;
    params
      .set<ParametersList>( "structureFunctions", structure_functions->parameters() )
      .set<int>( "mode", (int)mode )
      .set<int>( "beam1id", incoming_beams.first.pdg )
      .set<double>( "beam1pz", incoming_beams.first.pz )
      .set<int>( "beam2id", incoming_beams.second.pdg )
      .set<double>( "beam2pz", incoming_beams.second.pz )
      .set<std::vector<int> >( "ktFluxes", { (int)incoming_beams.first.kt_flux, (int)incoming_beams.second.kt_flux } )
      .set<double>( "sqrtS", sqrtS() );
    const HeavyIon hi1( incoming_beams.first.pdg ), hi2( incoming_beams.second.pdg );
    if ( hi1 )
      params.set<int>( "beam1A", hi1.A ).set<int>( "beam1Z", (int)hi1.Z );
    if ( hi2 )
      params.set<int>( "beam2A", hi2.A ).set<int>( "beam2Z", (int)hi2.Z );
    for ( auto& lim : cuts.central.list() ) {
      params.set<Limits>( lim.name, lim.limits );
      if ( lim.limits.hasMin() )
        params.set<double>( lim.name+"min", lim.limits.min() );
      if ( lim.limits.hasMax() )
        params.set<double>( lim.name+"max", lim.limits.max() );
    }
    if ( cuts.remnants.mass_single().valid() ) {
      params.set<Limits>( "mx", cuts.remnants.mass_single() );
      if ( cuts.remnants.mass_single().hasMin() )
        params.set<double>( "mxmin", cuts.remnants.mass_single().min() );
      if ( cuts.remnants.mass_single().hasMax() )
        params.set<double>( "mxmax", cuts.remnants.mass_single().max() );
    }
    if ( cuts.remnants.rapidity_single().valid() ) {
      params.set<Limits>( "yj", cuts.remnants.rapidity_single() );
      if ( cuts.remnants.rapidity_single().hasMin() )
        params.set<double>( "yjmin", cuts.remnants.rapidity_single().min() );
      if ( cuts.remnants.rapidity_single().hasMax() )
        params.set<double>( "yjmax", cuts.remnants.rapidity_single().max() );
    }
    if ( cuts.remnants.energy_single().valid() ) {
      const auto lim_xi = cuts.remnants.energy_single()*( -2./sqrtS() )+1.;
      params
        .set<Limits>( "xi", lim_xi )
        .set<double>( "ximin", lim_xi.hasMin() ? lim_xi.min() : 0. )
        .set<double>( "ximax", lim_xi.hasMax() ? lim_xi.max() : 1. );
    }

    return params;
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
