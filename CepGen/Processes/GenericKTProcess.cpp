#include "GenericKTProcess.h"
#include "CepGen/StructureFunctions/StructureFunctionsBuilder.h"
#include "CepGen/StructureFunctions/SigmaRatio.h"
#include "CepGen/Core/Exception.h"

namespace CepGen
{
  namespace Process
  {
    GenericKTProcess::GenericKTProcess( const std::string& name,
                                        const std::string& description,
                                        const std::array<ParticleCode,2>& partons,
                                        const std::vector<ParticleCode>& central ) :
      GenericProcess( name, description+" (kT-factorisation approach)" ),
      num_dimensions_( 0 ), kt_jacobian_( 0. ),
      qt1_( 0. ), phi_qt1_( 0. ),
      qt2_( 0. ), phi_qt2_( 0. ),
      flux1_( 0. ), flux2_( 0. ),
      kIntermediateParts( partons ), kProducedParts( central )
    {}

    void
    GenericKTProcess::addEventContent()
    {
      GenericProcess::setEventContent(
        { // incoming state
          { Particle::IncomingBeam1, Proton },
          { Particle::IncomingBeam2, Proton },
          { Particle::Parton1, kIntermediateParts[0] },
          { Particle::Parton2, kIntermediateParts[1] }
        },
        { // outgoing state
          { Particle::OutgoingBeam1, { Proton } },
          { Particle::OutgoingBeam2, { Proton } },
          { Particle::CentralSystem, kProducedParts }
        }
      );
      setExtraContent();
    }

    unsigned int
    GenericKTProcess::numDimensions( const Kinematics::ProcessMode& ) const
    {
      return num_dimensions_;
    }

    void
    GenericKTProcess::setKinematics( const Kinematics& kin )
    {
      cuts_ = kin;

      //============================================================================================
      // initialise the "constant" (wrt x) part of the Jacobian
      //============================================================================================

      kt_jacobian_ = 1.;
      num_dimensions_ = 0;
      mapped_variables_.clear();

      //============================================================================================
      // register the incoming partons' variables
      //============================================================================================

      registerVariable( qt1_, kLogarithmic, cuts_.cuts.initial[Cuts::qt], { 1.e-10, 500. }, "First incoming parton virtuality" );
      registerVariable( qt2_, kLogarithmic, cuts_.cuts.initial[Cuts::qt], { 1.e-10, 500. }, "Second incoming parton virtuality" );
      registerVariable( phi_qt1_, kLinear, cuts_.cuts.initial[Cuts::phi_qt], { 0., 2.*M_PI }, "First incoming parton azimuthal angle" );
      registerVariable( phi_qt2_, kLinear, cuts_.cuts.initial[Cuts::phi_qt], { 0., 2.*M_PI }, "Second incoming parton azimuthal angle" );

      //============================================================================================
      // register all process-dependent variables
      //============================================================================================

      preparePhaseSpace();

      //============================================================================================
      // register the outgoing remnants' variables
      //============================================================================================

      MX_ = MY_ = event_->getOneByRole( Particle::IncomingBeam1 ).mass();
      if ( cuts_.mode == Kinematics::InelasticElastic || cuts_.mode == Kinematics::InelasticInelastic )
        registerVariable( MX_, kSquare, cuts_.cuts.remnants[Cuts::mass], { 1.07, 1000. }, "Positive z proton remnant mass" );
      if ( cuts_.mode == Kinematics::ElasticInelastic || cuts_.mode == Kinematics::InelasticInelastic )
        registerVariable( MY_, kSquare, cuts_.cuts.remnants[Cuts::mass], { 1.07, 1000. }, "Negative z proton remnant mass" );
    }

    double
    GenericKTProcess::computeWeight()
    {
      if ( kt_jacobian_ == 0. )
        throw FatalError( "GenericKTProcess" )
          << "Point-independant component of the Jacobian for this "
          << "kt-factorised process is null.\n\tPlease check the "
          << "validity of the phase space!";

      //============================================================================================
      // generate and initialise all variables, and auxiliary (x-dependent) part of the Jacobian
      // for this phase space point.
      //============================================================================================

      const double aux_jacobian = generateVariables();
      if ( aux_jacobian == 0. )
        return 0.;

      //============================================================================================
      // compute the integrand and combine together into a single weight for the phase space point.
      //============================================================================================

      const double integrand = computeKTFactorisedMatrixElement();
      const double weight = ( kt_jacobian_*aux_jacobian ) * integrand;

      if ( Logger::get().level >= Logger::DebugInsideLoop )
        DebuggingInsideLoop( "GenericKTProcess" )
          << "Jacobian = " << kt_jacobian_ << " * " << aux_jacobian
          << "\n\tIntegrand = " << integrand
          << "\n\tdW = " << weight << ".";

      return weight;
    }

    const double GenericKTProcess::kMinFlux = 1.e-20;

    void
    GenericKTProcess::computeIncomingFluxes( double x1, double q1t2, double x2, double q2t2 )
    {
      flux1_ = flux2_ = 0.;
      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic:
          flux1_ = elasticFlux( x1, q1t2 );
          flux2_ = elasticFlux( x2, q2t2 );
          break;
        case Kinematics::ElasticInelastic:
          flux1_ = elasticFlux( x1, q1t2 );
          flux2_ = inelasticFlux( x2, q2t2, MY_, cuts_.structure_functions );
          break;
        case Kinematics::InelasticElastic:
          flux1_ = inelasticFlux( x1, q1t2, MX_, cuts_.structure_functions );
          flux2_ = elasticFlux( x2, q2t2 );
          break;
        case Kinematics::InelasticInelastic:
          flux1_ = inelasticFlux( x1, q1t2, MX_, cuts_.structure_functions );
          flux2_ = inelasticFlux( x2, q2t2, MY_, cuts_.structure_functions );
          break;
        default:
          throw FatalError( "GenericKTProcess" ) << "Invalid kinematics mode selected!";
      }
      flux1_ = std::max( flux1_, kMinFlux );
      flux2_ = std::max( flux2_, kMinFlux );
      DebuggingInsideLoop( "GenericKTProcess" ) << "Form factors: " << flux1_ << " / " << flux2_ << ".";
    }

    void
    GenericKTProcess::registerVariable( double& out, const MappingType& type,
                                        const Kinematics::Limits& in, Kinematics::Limits default_limits,
                                        const char* description )
    {
      Kinematics::Limits lim = in;
      out = 0.; // reset the variable
      if ( !in.valid() ) {
        Debugging( "GenericKTProcess" )
          << description << " could not be retrieved from the user configuration!\n\t"
          << "Setting it to the default value: " << default_limits << ".";
        lim = default_limits;
      }
      if ( type == kLogarithmic )
        lim = {
          std::max( log( lim.min() ), -10. ),
          std::min( log( lim.max() ), +10. )
        };
      mapped_variables_.emplace_back( MappingVariable{ lim, out, type, num_dimensions_++ } );
      switch ( type ) {
        case kSquare:
          kt_jacobian_ *= 2.*lim.range();
          break;
        default:
          kt_jacobian_ *= lim.range();
          break;
      }
      if ( Logger::get().level >= Logger::Debug )
        Debugging( "GenericKTProcess" )
          << description << " has been mapped to variable " << num_dimensions_ << ".\n\t"
          << "Allowed range for integration: " << lim << ".\n\t"
          << "Variable integration mode: " << type << ".";
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
          case kLinear: {
            cut.variable = cut.limits.x( xv );
          } break;
          case kLogarithmic: {
            cut.variable = exp( cut.limits.x( xv ) );
            jacobian *= cut.variable;
          } break;
          case kSquare: {
            cut.variable = cut.limits.x( xv );
            jacobian *= cut.variable;
          } break;
        }
      }
      if ( Logger::get().level >= Logger::DebugInsideLoop ) {
        std::ostringstream oss;
        for ( const auto& cut : mapped_variables_ ) {
          oss << "variable " << cut.index
              << " in range " << std::left << std::setw( 20 ) << cut.limits << std::right
              << " has value " << cut.variable << "\n\t";
        }
        DebuggingInsideLoop( "GenericKTProcess" ) << oss.str();
      }
      return jacobian;
    }

    void
    GenericKTProcess::fillKinematics( bool )
    {
      fillPrimaryParticlesKinematics();
      fillCentralParticlesKinematics(); // process-dependent!
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

      switch ( cuts_.mode ) {
        case Kinematics::ElasticElastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::ElasticInelastic:
          op1.setStatus( Particle::FinalState );
          op2.setStatus( Particle::Unfragmented ); op2.setMass( MY_ );
          break;
        case Kinematics::InelasticElastic:
          op1.setStatus( Particle::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::FinalState );
          break;
        case Kinematics::InelasticInelastic:
          op1.setStatus( Particle::Unfragmented ); op1.setMass( MX_ );
          op2.setStatus( Particle::Unfragmented ); op2.setMass( MY_ );
          break;
        default: {
          throw FatalError( "GenericKTProcess" )
            << "This kT factorisation process is intended for p-on-p collisions! Aborting.";
        } break;
      }

      //============================================================================================
      //     incoming partons (photons, pomerons, ...)
      //============================================================================================

      Particle& g1 = event_->getOneByRole( Particle::Parton1 );
      g1.setMomentum( event_->getOneByRole( Particle::IncomingBeam1 ).momentum()-PX_, true );
      g1.setStatus( Particle::Incoming );

      Particle& g2 = event_->getOneByRole( Particle::Parton2 );
      g2.setMomentum( event_->getOneByRole( Particle::IncomingBeam2 ).momentum()-PY_, true );
      g2.setStatus( Particle::Incoming );

      //============================================================================================
      //     two-parton system
      //============================================================================================

      event_->getOneByRole( Particle::Intermediate ).setMomentum( g1.momentum()+g2.momentum() );
    }

    double
    GenericKTProcess::elasticFlux( double x, double kt2 )
    {
      const double Q2_min = x*x*mp2_/( 1.-x ), Q2_ela = ( kt2+x*x*mp2_ )/( 1.-x );
      const FormFactors ff = FormFactors::ProtonElastic( Q2_ela );
      const double ela1 = ( 1.-x )*( 1.-Q2_min/Q2_ela );
      //const double ela3 = 1.-( Q2_ela-kt2 )/Q2_ela;

      const double f_ela = Constants::alphaEM*M_1_PI/kt2*( ela1*ff.FE + 0.5*x*x*ff.FM );

      // last factor below the Jacobian from dQ^2/Q^2 --> dkT^2/kT^2*(kT^2/Q^2)
      return f_ela*( 1.-x )*kt2 / Q2_ela;
    }

    double
    GenericKTProcess::inelasticFlux( double x, double kt2, double mx, const StructureFunctions::Type& sf, const FluxTypes& ft )
    {
      const double mx2 = mx*mx;

      // F2 structure function
      const double Q2min = 1. / ( 1.-x )*( x*( mx2-mp2_ ) + x*x*mp2_ ),
                   Q2 = Q2min + kt2/( 1.-x );
      float xbj = Q2 / ( Q2 + mx2 - mp2_ );

      const StructureFunctions str_fun = StructureFunctionsBuilder::get( sf, Q2, xbj );

      const double term1 = ( 1.-x )*( 1.-Q2min/Q2 );

      const double f_D = str_fun.F2/( mx2 + Q2 - mp2_ ) * term1;
      const double f_C = str_fun.F1( Q2, xbj ) * 2./Q2;

      return Constants::alphaEM*M_1_PI*( 1.-x )/Q2*( f_D+0.5*x*x*f_C );
    }

    std::ostream&
    operator<<( std::ostream& os, const GenericKTProcess::MappingType& type )
    {
      switch ( type ) {
        case GenericKTProcess::kLinear: return os << "linear";
        case GenericKTProcess::kLogarithmic: return os << "logarithmic";
        case GenericKTProcess::kSquare: return os << "squared";
      }
      return os;
    }
  }
}
