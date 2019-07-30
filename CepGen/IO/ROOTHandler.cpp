#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Parameters.h"

#include "CepGen/Version.h"

#include "TFile.h"
#include "TH1.h"

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the generic ROOT file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class ROOTHandler : public GenericExportHandler
    {
      public:
        explicit ROOTHandler( const ParametersList& );
        ~ROOTHandler();

        void initialise( const Parameters& ) override {}
        void setCrossSection( double xsec, double ) override { xsec_ = xsec; }
        void operator<<( const Event& ) override;

      private:
        std::unique_ptr<TFile> file_;
        std::vector<std::pair<std::string,TH1*> > hists_;
        const ParametersList variables_;

        double xsec_;
        const utils::EventBrowser browser_;
    };

    ROOTHandler::ROOTHandler( const ParametersList& params ) :
      GenericExportHandler( "root" ),
      file_     ( TFile::Open( params.get<std::string>( "filename", "output.root" ).c_str(), "recreate" ) ),
      variables_( params.get<ParametersList>( "variables" ) ),
      xsec_( 1. )
    {
      //--- extract list of variables to be plotted in histograms
      for ( const auto& var : variables_.keys() ) {
        const auto& hvar = variables_.get<ParametersList>( var );
        const int nbins = hvar.get<int>( "nbins", 10 );
        const double min = hvar.get<double>( "low", 0. ), max = hvar.get<double>( "high", 1. );
        const auto title = Form( "%s;%s;d#sigma/d(%s) (pb/bin)", var.c_str(), var.c_str(), var.c_str() );
        hists_.emplace_back( std::make_pair( var, new TH1D( var.c_str(), title.c_str(), nbins, min, max ) ) );
        CG_INFO( "ROOTHandler" )
          << "Booking a histogram with " << nbins << " bin" << utils::s( nbins )
          << " between " << min << " and " << max << " for \"" << var << "\".";
      }
    }

    ROOTHandler::~ROOTHandler()
    {
      //--- finalisation of the output file
      for ( const auto& hist : hists_ )
        hist.second->Write( hist.first.c_str() );
        // ROOT and its sumptuous memory management disallows the "delete" here
      file_->Close();
    }

    void
    ROOTHandler::operator<<( const Event& ev )
    {
      //--- increment the corresponding histograms
      for ( const auto& h_var : hists_ )
        h_var.second->Fill( browser_.get( ev, h_var.first ), xsec_ );
    }
  }
}

REGISTER_IO_MODULE( "root", ROOTHandler )
