#include "CepGen/Processes/GenericKTProcess.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"

#include "CepGen/Physics/Constants.h"
#include "CepGen/Physics/KTFlux.h"
#include "CepGen/Physics/FormFactors.h"
#include "CepGen/Physics/PDG.h"

#include <iomanip>

namespace cepgen
{
  namespace proc
  {
    GenericKTProcess::GenericKTProcess( const ParametersList& params,
                                        const std::string& name,
                                        const std::string& description,
                                        const std::array<pdgid_t,2>& partons,
                                        const std::vector<pdgid_t>& central ) :
      GenericProcess( params, name, description+" (kT-factorisation approach)" ),
      num_dimensions_( 0 ), kt_jacobian_( 0. ),
      qt1_( 0. ), phi_qt1_( 0. ), qt2_( 0. ), phi_qt2_( 0. ),
      kIntermediateParts( partons ), kProducedParts( central )
    {}

    void
    GenericKTProcess::addEventContent()
    {
      GenericProcess::setEventContent(
        { // incoming state
          { Particle::IncomingBeam1, PDG::proton },
          { Particle::IncomingBeam2, PDG::proton },
          { Particle::Parton1, kIntermediateParts[0] },
          { Particle::Parton2, kIntermediateParts[1] }
        },
        { // outgoing state
          { Particle::OutgoingBeam1, { PDG::proton } },
          { Particle::OutgoingBeam2, { PDG::proton } },
          { Particle::CentralSystem, kProducedParts }
        }
      );
      setExtraContent();
    }

    unsigned int
    GenericKTProcess::numDimensions() const
    {
      return num_dimensions_;
    }

    void
    GenericKTProcess::setKinematics( const Kinematics& kin )
    {
      kin_ = kin;

      const KTFlux flux1 = (KTFlux)kin_.incoming_beams.first.kt_flux,
                   flux2 = (KTFlux)kin_.incoming_beams.second.kt_flux;

      if ( kin_.mode == KinematicsMode::invalid ) {
        //--- try to extrapolate kinematics mode from unintegrated fluxes
        bool el1 = ( flux1 == KTFlux::P_Photon_Elastic
                  || flux1 == KTFlux::HI_Photon_Elastic
                  || flux1 == KTFlux::P_Gluon_KMR );
        bool el2 = ( flux2 == KTFlux::P_Photon_Elastic
                  || flux2 == KTFlux::HI_Photon_Elastic
                  || flux2 == KTFlux::P_Gluon_KMR );
        if ( el1 && el2 )
          kin_.mode = KinematicsMode::ElasticElastic;
        else if ( el1 )
          kin_.mode = KinematicsMode::ElasticInelastic;
        else if ( el2 )
          kin_.mode = KinematicsMode::InelasticElastic;
        else
          kin_.mode = KinematicsMode::InelasticInelastic;
      }
      else {
        //--- try to extrapolate unintegrated fluxes from kinematics mode
        const HeavyIon hi1( kin_.incoming_beams.first.pdg ), hi2( kin_.incoming_beams.second.pdg );
        //==========================================================================================
        // ensure the first incoming flux is compatible with the kinematics mode
        //==========================================================================================
        if ( ( kin_.mode == KinematicsMode::ElasticElastic
            || kin_.mode == KinematicsMode::ElasticInelastic )
          && ( flux1 != KTFlux::P_Photon_Elastic ) ) {
          kin_.incoming_beams.first.kt_flux = hi1
            ? KTFlux::HI_Photon_Elastic
            : KTFlux::P_Photon_Elastic;
          CG_DEBUG( "GenericKTProcess:kinematics" )
            << "Set the kt flux for first incoming photon to \""
            << kin_.incoming_beams.first.kt_flux << "\".";
        }
        else if ( flux1 != KTFlux::P_Photon_Inelastic
               && flux1 != KTFlux::P_Photon_Inelastic_Budnev ) {
          if ( hi1 )
            throw CG_FATAL( "GenericKTProcess:kinematics" )
              << "Inelastic photon emission from HI not yet supported!";
          kin_.incoming_beams.first.kt_flux = KTFlux::P_Photon_Inelastic_Budnev;
          CG_DEBUG( "GenericKTProcess:kinematics" )
            << "Set the kt flux for first incoming photon to \""
            << kin_.incoming_beams.first.kt_flux << "\".";
        }
        //==========================================================================================
        // ensure the second incoming flux is compatible with the kinematics mode
        //==========================================================================================
        if ( ( kin_.mode == KinematicsMode::ElasticElastic
            || kin_.mode == KinematicsMode::InelasticElastic )
          && ( flux2 != KTFlux::P_Photon_Elastic ) ) {
          kin_.incoming_beams.second.kt_flux = hi2
            ? KTFlux::HI_Photon_Elastic
            : KTFlux::P_Photon_Elastic;
          CG_DEBUG( "GenericKTProcess:kinematics" )
            << "Set the kt flux for second incoming photon to \""
            << kin_.incoming_beams.second.kt_flux << "\".";
        }
        else if ( flux2 != KTFlux::P_Photon_Inelastic
               && flux2 != KTFlux::P_Photon_Inelastic_Budnev ) {
          if ( hi2 )
            throw CG_FATAL( "GenericKTProcess:kinematics" )
              << "Inelastic photon emission from HI not yet supported!";
          kin_.incoming_beams.second.kt_flux = KTFlux::P_Photon_Inelastic_Budnev;
          CG_DEBUG( "GenericKTProcess:kinematics" )
            << "Set the kt flux for second incoming photon to \""
            << kin_.incoming_beams.second.kt_flux << "\".";
        }
      }

      //============================================================================================
      // initialise the "constant" (wrt x) part of the Jacobian
      //============================================================================================

      kt_jacobian_ = 1.;
      num_dimensions_ = 0;
      mapped_variables_.clear();

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      registerVariable( qt1_, Mapping::logarithmic, kin_.cuts.initial.qt, { 1.e-10, 500. }, "First incoming parton virtuality" );
      registerVariable( qt2_, Mapping::logarithmic, kin_.cuts.initial.qt, { 1.e-10, 500. }, "Second incoming parton virtuality" );
      registerVariable( phi_qt1_, Mapping::linear, kin_.cuts.initial.phi_qt, { 0., 2.*M_PI }, "First incoming parton azimuthal angle" );
      registerVariable( phi_qt2_, Mapping::linear, kin_.cuts.initial.phi_qt, { 0., 2.*M_PI }, "Second incoming parton azimuthal angle" );

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      MX_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
      MY_ = event_->getOneByRole( Particle::IncomingBeam2 ).mass();
      if ( kin_.mode == KinematicsMode::InelasticElastic || kin_.mode == KinematicsMode::InelasticInelastic )
        registerVariable( MX_, Mapping::square, kin_.cuts.remnants.mass_single, { 1.07, 1000. }, "Positive z proton remnant mass" );
      if ( kin_.mode == KinematicsMode::ElasticInelastic || kin_.mode == KinematicsMode::InelasticInelastic )
        registerVariable( MY_, Mapping::square, kin_.cuts.remnants.mass_single, { 1.07, 1000. }, "Negative z proton remnant mass" );

      prepareKinematics();
    }

    double
    GenericKTProcess::computeWeight()
    {
      if ( mapped_variables_.size() == 0 )
        throw CG_FATAL( "GenericKTProcess:weight" )
          << "No variables are mapped with this process!";
      if ( kt_jacobian_ == 0. )
        throw CG_FATAL( "GenericKTProcess:weight" )
          << "Point-independant component of the Jacobian for this "
          << "kt-factorised process is null.\n\t"
          << "Please check the validity of the phase space!";

      //============================================================================================
      // generate and initialise all variables, and auxiliary (x-dependent) part of the Jacobian
      // for this phase space point.
      //============================================================================================

      const double aux_jacobian = generateVariables();
      if ( aux_jacobian <= 0. )
        return 0.;

      //============================================================================================
      // compute the integrand and combine together into a single weight for the phase space point.
      //============================================================================================

      const double integrand = computeKTFactorisedMatrixElement();
      if ( integrand <= 0. )
        return 0.;

      const double weight = ( kt_jacobian_*aux_jacobian ) * integrand;

      CG_DEBUG_LOOP( "GenericKTProcess:weight" )
        << "Jacobian: " << kt_jacobian_ << " * " << aux_jacobian
        << " = " << ( kt_jacobian_*aux_jacobian ) << ".\n\t"
        << "Integrand = " << integrand << "\n\t"
        << "dW = " << weight << ".";

      return weight;
    }

    void
    GenericKTProcess::registerVariable( double& out, const Mapping& type,
                                        const Limits& in, Limits default_limits,
                                        const char* description )
    {
      Limits lim = in;
      out = 0.; // reset the variable
      if ( !in.valid() ) {
        CG_DEBUG( "GenericKTProcess:registerVariable" )
          << description << " could not be retrieved from the user configuration!\n\t"
          << "Setting it to the default value: " << default_limits << ".";
        lim = default_limits;
      }
      if ( type == Mapping::logarithmic )
        lim = {
          std::max( log( lim.min() ), -10. ),
          std::min( log( lim.max() ), +10. )
        };
      mapped_variables_.emplace_back( MappingVariable{ description, lim, out, type, num_dimensions_++ } );
      switch ( type ) {
        case Mapping::square:
          kt_jacobian_ *= 2.*lim.range();
          break;
        default:
          kt_jacobian_ *= lim.range();
          break;
      }
      CG_DEBUG( "GenericKTProcess:registerVariable" )
        << description << " has been mapped to variable " << num_dimensions_ << ".\n\t"
        << "Allowed range for integration: " << lim << ".\n\t"
        << "Variable integration mode: " << type << ".";
    }

    void
    GenericKTProcess::dumpVariables() const
    {
      std::ostringstream os;
      for ( const auto& var : mapped_variables_ )
        os << "\n\t(" << var.index << ") " << var.type << " mapping (" << var.description << ") in range " << var.limits;
      CG_INFO( "GenericKTProcess:dumpVariables" )
        << "List of variables handled by this kt-factorised process:"
        << os.str();
    }

    double
    GenericKTProcess::generateVariables() const
    {
      double jacobian = 1.;
      for ( const auto& cut : mapped_variables_ ) {
        if ( !cut.limits.valid() )
          continue;
        const double xv = x( cut.index ); // between 0 and 1
        switch ( cut.type ) {
          case Mapping::linear: {
            cut.variable = cut.limits.x( xv );
          } break;
          case Mapping::logarithmic: {
            cut.variable = exp( cut.limits.x( xv ) );
            jacobian *= cut.variable;
          } break;
          case Mapping::square: {
            cut.variable = cut.limits.x( xv );
            jacobian *= cut.variable;
          } break;
        }
      }
      if ( CG_LOG_MATCH( "KtProcess:vars", debugInsideLoop ) ) {
        std::ostringstream oss;
        for ( const auto& cut : mapped_variables_ ) {
          oss << "variable " << cut.index
              << " in range " << std::left << std::setw( 20 ) << cut.limits << std::right
              << " has value " << cut.variable << "\n\t";
        }
        CG_DEBUG_LOOP( "KtProcess:vars" ) << oss.str();
      }
      return jacobian;
    }

    void
    GenericKTProcess::fillKinematics( bool )
    {
      fillCentralParticlesKinematics(); // process-dependent!
      fillPrimaryParticlesKinematics();
    }

    void
    GenericKTProcess::fillPrimaryParticlesKinematics()
    {
      //============================================================================================
      //     outgoing protons
      //============================================================================================

      Particle& op1 = event_->getOneByRole( Particle::OutgoingBeam1 ),
               &op2 = event_->getOneByRole( Particle::OutgoingBeam2 );

      op1.setMomentum( PX_ );
      op2.setMomentum( PY_ );

      switch ( kin_.mode ) {
        case KinematicsMode::ElasticElastic:
          op1.setStatus( Particle::Status::FinalState );
          op2.setStatus( Particle::Status::FinalState );
          break;
        case KinematicsMode::ElasticInelastic:
          op1.setStatus( Particle::Status::FinalState );
          op2.setStatus( Particle::Status::Unfragmented ); op2.setMass( MY_ );
          break;
        case KinematicsMode::InelasticElastic:
          op1.setStatus( Particle::Status::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::Status::FinalState );
          break;
        case KinematicsMode::InelasticInelastic:
          op1.setStatus( Particle::Status::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::Status::Unfragmented ); op2.setMass( MY_ );
          break;
        default: {
          throw CG_FATAL( "GenericKTProcess" )
            << "This kT factorisation process is intended for p-on-p collisions! Aborting.";
        } break;
      }

      //============================================================================================
      //     incoming partons (photons, pomerons, ...)
      //============================================================================================

      Particle& g1 = event_->getOneByRole( Particle::Parton1 );
      g1.setMomentum( event_->getOneByRole( Particle::IncomingBeam1 ).momentum()-PX_, true );

      Particle& g2 = event_->getOneByRole( Particle::Parton2 );
      g2.setMomentum( event_->getOneByRole( Particle::IncomingBeam2 ).momentum()-PY_, true );

      //============================================================================================
      //     two-parton system
      //============================================================================================

      event_->getOneByRole( Particle::Intermediate ).setMomentum( g1.momentum()+g2.momentum() );
    }

    std::ostream&
    operator<<( std::ostream& os, const GenericKTProcess::Mapping& type )
    {
      switch ( type ) {
        case GenericKTProcess::Mapping::linear: return os << "linear";
        case GenericKTProcess::Mapping::logarithmic: return os << "logarithmic";
        case GenericKTProcess::Mapping::square: return os << "squared";
      }
      return os;
    }
  }
}
