#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Parameters.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"

#include "modules/Delphes.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "TFile.h"

#include <sstream>

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Export handler for Delphes
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class DelphesHandler : public GenericExportHandler
    {
      public:
        struct CepGenConfReader : public ExRootConfReader
        {
          using ExRootConfReader::ExRootConfReader;
          void feedParameters( const Parameters& );
        };
        /// Class constructor
        explicit DelphesHandler( const ParametersList& );
        ~DelphesHandler();

        void initialise( const Parameters& ) override;
        /// Writer operator
        void operator<<( const Event& ) override;
        void setCrossSection( double, double ) override {}

      private:
        std::unique_ptr<TFile> output_;
        const std::string input_card_;
        std::unique_ptr<Delphes> delphes_;
        std::unique_ptr<ExRootConfReader> conf_reader_;
        std::unique_ptr<ExRootTreeWriter> tree_writer_;
        TObjArray* out_all_parts_, *out_stab_parts_, *out_partons_;
        DelphesFactory* factory_;
    };

    DelphesHandler::DelphesHandler( const ParametersList& params ) :
      GenericExportHandler( "delphes" ),
      output_( new TFile( params.get<std::string>( "filename", "output.delphes.root" ).c_str(), "recreate" ) ),
      input_card_( params.get<std::string>( "inputCard", "input.tcl" ) ),
      delphes_( new Delphes ),
      conf_reader_( new ExRootConfReader ),
      tree_writer_( new ExRootTreeWriter( output_.get(), "Delphes") ),
      out_all_parts_( nullptr ), out_stab_parts_( nullptr ), out_partons_( nullptr ),
      factory_( nullptr )
    {
      conf_reader_->ReadFile( input_card_.c_str() );
      delphes_->SetTreeWriter( tree_writer_.get() );
      delphes_->SetConfReader( conf_reader_.get() );
    }

    DelphesHandler::~DelphesHandler()
    {
      delphes_->FinishTask();
      tree_writer_->Write();
    }

    void
    DelphesHandler::initialise( const Parameters& params )
    {
      /*CepGenConfReader conf;
      conf.feedParameters( params );
      delphes_->SetConfReader( &conf );*/
      factory_ = delphes_->GetFactory();
      out_all_parts_ = delphes_->ExportArray( "allParticles" );
      out_stab_parts_ = delphes_->ExportArray( "stableParticles" );
      out_partons_ = delphes_->ExportArray( "partons" );
      delphes_->InitTask();
    }

    void
    DelphesHandler::operator<<( const Event& ev )
    {
      delphes_->Clear();
      tree_writer_->Clear();
      //...
      for ( const auto& part : ev.particles() ) {
        auto cand = factory_->NewCandidate();
        cand->PID = part.integerPdgId();
        cand->Status = (int)part.status();
        cand->Charge = part.charge();
        cand->Mass = part.mass();
        const auto& mom = part.momentum();
        cand->Momentum.SetPxPyPzE( mom.px(), mom.py(), mom.pz(), mom.energy() );
        cand->M1 = part.primary() ? 0 : *part.mothers().begin();
        cand->M2 = part.mothers().size() < 2 ? 0 : *part.mothers().rbegin();
        cand->D1 = part.daughters().empty() ? -1 : *part.daughters().begin();
        cand->D2 = part.daughters().size() < 2 ? -1 : *part.daughters().rbegin();
        out_all_parts_->Add( cand );
        if ( cand->Status == 1 )
          out_stab_parts_->Add( cand );
      }
      //...
      delphes_->ProcessTask();
      tree_writer_->Fill();
    }

    void
    DelphesHandler::CepGenConfReader::feedParameters( const Parameters& params )
    {
    }
  }
}

REGISTER_IO_MODULE( delphes, DelphesHandler )
