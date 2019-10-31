#include "CepGen/Core/GenericExportHandler.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/EventModifier.h"
#include "CepGen/Core/utils.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include <sstream>

namespace cepgen
{
  namespace io
  {
    GenericExportHandler::GenericExportHandler( const ParametersList& params ) :
      params_( params ),
      name_( params_.get<std::string>( "mod_name" ) ),
      event_num_( 0. )
    {}

    GenericExportHandler::~GenericExportHandler()
    {}

    std::string
    GenericExportHandler::banner( const Parameters& params, const std::string& prep )
    {
      std::ostringstream os;
      os
        << prep << "  ***** Sample generated with CepGen v" << version() << " *****\n"
        << prep << "  * process: " << params.processName() << " (" << params.kinematics.mode << ")\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic ) {
        os << prep << "  * structure functions: " << params.kinematics.structure_functions->description() << "\n";
        if ( !params.eventModifiersSequence().empty() ) {
          os << prep << "  * " << utils::s( "event modifier", params.eventModifiersSequence().size() ) << ": ";
          std::string sep;
          for ( const auto& mod : params.eventModifiersSequence() )
            os << sep << mod->name(), sep = ", ";
          os << "\n";
        }
      }
      os
        << prep << "  *--- incoming state\n";
      if ( params.kinematics.cuts.initial.q2.valid() )
        os
          << prep << "  * Q2 range (GeV2): "
          << params.kinematics.cuts.initial.q2 << "\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic
        && params.kinematics.cuts.remnants.mass_single.valid() )
        os
          << prep << "  * remnants mass range (GeV/c2): "
          << params.kinematics.cuts.remnants.mass_single << "\n";
      os << prep << "  *--- central system\n";
      if ( params.kinematics.cuts.central.pt_single.valid() )
        os
          << prep << "  * single particle pt (GeV/c): "
          << params.kinematics.cuts.central.pt_single << "\n";
      if ( params.kinematics.cuts.central.energy_single.valid() )
        os
          << prep << "  * single particle energy (GeV): "
          << params.kinematics.cuts.central.energy_single << "\n";
      if ( params.kinematics.cuts.central.eta_single.valid() )
        os
          << prep << "  * single particle eta: "
          << params.kinematics.cuts.central.eta_single << "\n";
      if ( params.kinematics.cuts.central.pt_sum.valid() )
        os
          << prep << "  * total pt (GeV/c): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      if ( params.kinematics.cuts.central.mass_sum.valid() )
        os
          << prep << "  * total invariant mass (GeV/c2): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      os
        << prep << "  **************************************************";
      return os.str();
    }
  }
}

