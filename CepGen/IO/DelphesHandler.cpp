#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Parameters.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"

#include "modules/Delphes.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesClasses.h"

#include <sstream>

namespace cepgen
{
  namespace output
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
        std::unique_ptr<Delphes> delphes_;
    };

    DelphesHandler::DelphesHandler( const ParametersList& params ) :
      delphes_( new Delphes )
    {}

    DelphesHandler::~DelphesHandler()
    {
      delphes_->FinishTask();
    }

    void
    DelphesHandler::initialise( const Parameters& params )
    {
      CepGenConfReader conf;
      conf.feedParameters( params );
      delphes_->SetConfReader( &conf );
      delphes_->InitTask();
    }

    void
    DelphesHandler::operator<<( const Event& ev )
    {
      delphes_->Clear();
      auto factory = delphes_->GetFactory();
      //...
      for ( const auto& part : ev.particles() ) {
        auto cand = factory->NewCandidate();
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
      }
      //...
      delphes_->ProcessTask();
    }

    void
    DelphesHandler::CepGenConfReader::feedParameters( const Parameters& params )
    {
    }
  }
}

REGISTER_IO_MODULE( delphes, DelphesHandler )
