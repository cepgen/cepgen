#include "CepGenAddOns/EventInterfaces/ROOTTreeInfo.h"

#include "CepGen/Modules/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Parameters.h"

// ROOT includes
#include "TFile.h"

#include <sstream>

namespace cepgen
{
  namespace io
  {
    /**
     * Handler for the storage of events in a ROOT format
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date 27 Jan 2014
     */
    class ROOTTreeHandler : public ExportModule
    {
      public:
        /// Class constructor
        explicit ROOTTreeHandler( const ParametersList& );
        ~ROOTTreeHandler();

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override;

      private:
        std::unique_ptr<TFile> file_;
        const bool compress_;
        std::unique_ptr<ROOT::CepGenRun> run_tree_;
        std::unique_ptr<ROOT::CepGenEvent> evt_tree_;
    };

    ROOTTreeHandler::ROOTTreeHandler( const ParametersList& params ) :
      ExportModule( params ),
      file_( TFile::Open( params.get<std::string>( "filename", "output.root" ).c_str(), "recreate" ) ),
      compress_( params.get<bool>( "compress", false ) ),
      run_tree_( new ROOT::CepGenRun ), evt_tree_( new ROOT::CepGenEvent )
    {
      if ( !file_->IsOpen() )
        throw CG_FATAL( "ROOTTreeHandler" ) << "Failed to create the output file!";
      run_tree_->create();
      evt_tree_->create();
    }

    ROOTTreeHandler::~ROOTTreeHandler()
    {
      run_tree_->fill();
      file_->Write();
    }

    void
    ROOTTreeHandler::initialise( const Parameters& params )
    {
      run_tree_->litigious_events = 0;
      run_tree_->sqrt_s = params.kinematics.sqrtS();
    }

    void
    ROOTTreeHandler::operator<<( const Event& ev )
    {
      evt_tree_->fill( ev, compress_ );
      run_tree_->num_events += 1;
    }

    void
    ROOTTreeHandler::setCrossSection( double xsect, double xsect_err )
    {
      run_tree_->xsect = xsect;
      run_tree_->errxsect = xsect_err;
    }
  }
}

REGISTER_IO_MODULE( "root_tree", ROOTTreeHandler )
