#include "CepGen/IO/ExportHandler.h"
#include "CepGen/IO/ROOTTreeInfo.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"

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
    class ROOTTreeHandler : public GenericExportHandler
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
        std::unique_ptr<ROOT::CepGenRun> run_tree_;
        std::unique_ptr<ROOT::CepGenEvent> evt_tree_;
    };

    ROOTTreeHandler::ROOTTreeHandler( const ParametersList& params ) :
      GenericExportHandler( "root" ),
      file_( TFile::Open( params.get<std::string>( "filename", "output.root" ).c_str(), "recreate" ) ),
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
      evt_tree_->gen_time = ev.time_generation;
      evt_tree_->tot_time = ev.time_total;
      evt_tree_->np = 0;
      for ( const auto& p : ev.particles() ) {
        const auto& m = p.momentum();
        evt_tree_->rapidity[evt_tree_->np] = m.rapidity();
        evt_tree_->pt[evt_tree_->np] = m.pt();
        evt_tree_->eta[evt_tree_->np] = m.eta();
        evt_tree_->phi[evt_tree_->np] = m.phi();
        evt_tree_->E[evt_tree_->np] = p.energy();
        evt_tree_->m[evt_tree_->np] = p.mass();
        evt_tree_->pdg_id[evt_tree_->np] = p.integerPdgId();
        evt_tree_->parent1[evt_tree_->np] = ( p.mothers().size() > 0 ) ? *p.mothers().begin() : -1;
        evt_tree_->parent2[evt_tree_->np] = ( p.mothers().size() > 1 ) ? *p.mothers().rbegin() : -1;
        evt_tree_->status[evt_tree_->np] = (int)p.status();
        evt_tree_->stable[evt_tree_->np] = ( (short)p.status() > 0 );
        evt_tree_->charge[evt_tree_->np] = p.charge();
        evt_tree_->role[evt_tree_->np] = p.role();

        evt_tree_->np++;
      }
      run_tree_->num_events += 1;
      evt_tree_->fill();
    }

    void
    ROOTTreeHandler::setCrossSection( double xsect, double xsect_err )
    {
      run_tree_->xsect = xsect;
      run_tree_->errxsect = xsect_err;
    }
  }
}

REGISTER_IO_MODULE( root_tree, ROOTTreeHandler )
