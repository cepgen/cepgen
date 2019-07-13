#include "CepGen/IO/GenericExportHandler.h"

#include "CepGen/StructureFunctions/StructureFunctions.h"
#include "CepGen/Physics/Constants.h"

#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include <sstream>

namespace cepgen
{
  namespace io
  {
    GenericExportHandler::GenericExportHandler( const std::string& name ) :
      name_( name ), event_num_( 0. )
    {}

    GenericExportHandler::~GenericExportHandler()
    {}

    std::string
    GenericExportHandler::banner( const Parameters& params )
    {
      std::ostringstream os;
      os
        << "  ***** Sample generated with CepGen v" << version() << " *****\n"
        << "  * process: " << params.processName() << " (" << params.kinematics.mode << ")\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic ) {
        os << "  * structure functions: " << params.kinematics.structure_functions->description() << "\n";
        if ( !params.hadroniserName().empty() )
          os << "  * hadroniser: " << params.hadroniserName() << "\n";
      }
      os
        << "  *--- incoming state\n";
      if ( params.kinematics.cuts.initial.q2.valid() )
        os
          << "  * Q2 range (GeV2): "
          << params.kinematics.cuts.initial.q2 << "\n";
      if ( params.kinematics.mode != KinematicsMode::ElasticElastic
        && params.kinematics.cuts.remnants.mass_single.valid() )
        os
          << "  * remnants mass range (GeV/c2): "
          << params.kinematics.cuts.remnants.mass_single << "\n";
      os << "  *--- central system\n";
      if ( params.kinematics.cuts.central.pt_single.valid() )
        os
          << "  * single particle pt (GeV/c): "
          << params.kinematics.cuts.central.pt_single << "\n";
      if ( params.kinematics.cuts.central.energy_single.valid() )
        os
          << "  * single particle energy (GeV): "
          << params.kinematics.cuts.central.energy_single << "\n";
      if ( params.kinematics.cuts.central.eta_single.valid() )
        os
          << "  * single particle eta: "
          << params.kinematics.cuts.central.eta_single << "\n";
      if ( params.kinematics.cuts.central.pt_sum.valid() )
        os
          << "  * total pt (GeV/c): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      if ( params.kinematics.cuts.central.mass_sum.valid() )
        os
          << "  * total invariant mass (GeV/c2): "
          << params.kinematics.cuts.central.mass_sum << "\n";
      os
        << "  **************************************************";
      return os.str();
    }
  }
}

