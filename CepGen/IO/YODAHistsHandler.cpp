#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"

#include "YODA/Histo1D.h"
#include "YODA/WriterYODA.h"
#include "YODA/WriterAIDA.h"
#include "YODA/WriterFLAT.h"

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the generic YODA file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class YODAHistsHandler : public GenericExportHandler
    {
      public:
        explicit YODAHistsHandler( const ParametersList& );
        ~YODAHistsHandler();

        void initialise( const Parameters& ) override {}
        void setCrossSection( double xsec, double ) override { xsec_ = xsec; }
        void operator<<( const Event& ) override;

      private:
        std::ofstream file_;
        std::vector<std::pair<std::string,YODA::Histo1D> > hists_;
        const std::string type_;
        const ParametersList variables_;

        double xsec_;
        const utils::EventBrowser browser_;
    };

    YODAHistsHandler::YODAHistsHandler( const ParametersList& params ) :
      GenericExportHandler( "yoda" ),
      file_( params.get<std::string>( "filename", "output.yoda" ) ),
      type_( params.get<std::string>( "type", "yoda" ) ),
      variables_( params.get<ParametersList>( "variables" ) ),
      xsec_( 1. )
    {
      //--- extract list of variables to be plotted in histograms
      for ( const auto& var : variables_.keys() ) {
        const auto& hvar = variables_.get<ParametersList>( var );
        const int nbins = hvar.get<int>( "nbins", 10 );
        const double min = hvar.get<double>( "low", 0. ), max = hvar.get<double>( "high", 1. );
        const auto title = Form( "d(sigma)/d(%s) (pb/bin)", var.c_str() );
        hists_.emplace_back( std::make_pair( var, YODA::Histo1D( nbins, min, max, var, title ) ) );
        CG_DEBUG( "YODAHistsHandler" )
          << "Booking a histogram with " << nbins << " bin" << utils::s( nbins )
          << " between " << min << " and " << max << " for \"" << var << "\".";
      }
    }

    YODAHistsHandler::~YODAHistsHandler()
    {
      //--- finalisation of the output file
      for ( const auto& hist : hists_ )
        if ( type_ == "yoda" )
          YODA::WriterYODA::create().write( file_, hist.second );
        else if ( type_ == "aida" )
          YODA::WriterAIDA::create().write( file_, hist.second );
        else if ( type_ == "flat" )
          YODA::WriterFLAT::create().write( file_, hist.second );
    }

    void
    YODAHistsHandler::operator<<( const Event& ev )
    {
      //--- increment the corresponding histograms
      for ( auto& h_var : hists_ )
        h_var.second.fillBin( browser_.get( ev, h_var.first ), xsec_ );
    }
  }
}

REGISTER_IO_MODULE( "yoda", YODAHistsHandler )
