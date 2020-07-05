#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Utils/String.h"

#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Counter.h"

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
     * \tparam T YODA writer type
     */
    template<typename T>
    class YODAHistsHandler : public ExportModule
    {
      public:
        explicit YODAHistsHandler( const ParametersList& );
        ~YODAHistsHandler();
        static std::string description() { return "YODA histograms/profiles file output module"; }

        void initialise( const Parameters& ) override {}
        void setCrossSection( double xsec, double ) override { xsec_ = xsec; }
        void operator<<( const Event& ) override;

      private:
        std::ofstream file_;
        std::vector<std::pair<std::string,YODA::Histo1D> > hists1d_;
        std::vector<std::pair<std::vector<std::string>,YODA::Histo2D> > hists2d_;
        std::vector<std::pair<std::vector<std::string>,YODA::Profile1D> > profiles1d_;
        std::vector<std::pair<std::vector<std::string>,YODA::Profile2D> > profiles2d_;
        YODA::Counter weight_cnt_;
        const ParametersList variables_;

        double xsec_;
        const utils::EventBrowser browser_;
    };

    template<typename T>
    YODAHistsHandler<T>::YODAHistsHandler( const ParametersList& params ) :
      ExportModule( params ),
      file_( params.get<std::string>( "filename", "output.yoda" ) ),
      variables_( params.get<ParametersList>( "variables" ) ),
      xsec_( 1. )
    {
      //--- extract list of variables/correlations to be plotted in histograms
      for ( const auto& key : variables_.keys() ) {
        const auto& vars = utils::split( key, ':' );
        if ( vars.size() < 1 || vars.size() > 3 )
          throw CG_FATAL( "YODAHistsHandler" )
            << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvars = variables_.get<ParametersList>( key );
        int nbins_x = hvars.get<int>( "nbinsX", 10 );
        nbins_x = hvars.get<int>( "nbins", nbins_x );
        double min_x = hvars.get<double>( "lowX", 0. ), max_x = hvars.get<double>( "highX", 1. );
        min_x = hvars.get<double>( "low", min_x ), max_x = hvars.get<double>( "high", max_x );
        const bool profile = hvars.get<bool>( "profile", false );
        if ( vars.size() == 1 ) { // 1D histogram
          const auto title = utils::format( "d(sigma)/d(%s) (pb/bin)", key.c_str() );
          hists1d_.emplace_back( std::make_pair( key, YODA::Histo1D( nbins_x, min_x, max_x, key, title ) ) );
          CG_INFO( "YODAHistsHandler" )
            << "Booking a histogram with " << utils::s( "bin", nbins_x )
            << " between " << min_x << " and " << max_x << " for \"" << vars[0] << "\".";
          continue;
        }
        const int nbins_y = hvars.get<int>( "nbinsY", 10 );
        const double min_y = hvars.get<double>( "lowY", 0. ), max_y = hvars.get<double>( "highY", 1. );
        if ( vars.size() == 2 ) { // 2D histogram
          const auto title = utils::format( "d^2(sigma)/d(%s)/d(%s) (pb/bin)", vars[0].c_str(), vars[1].c_str() );
          if ( profile )
            profiles1d_.emplace_back( std::make_pair( vars, YODA::Profile1D( nbins_x, min_x, max_x, key, title ) ) );
          else
            hists2d_.emplace_back( std::make_pair( vars, YODA::Histo2D( nbins_x, min_x, max_x, nbins_y, min_y, max_y, key, title ) ) );
          CG_INFO( "YODAHistsHandler" )
            << "Booking a "
            << ( profile ? "1D profile" : "2D correlation plot" )
            << " with " << utils::s( "bin", nbins_x+nbins_y )
            << " between (" << min_x << ", " << min_y << ") and (" << max_x << ", " << max_y << ")"
            << " for \"" << utils::merge( vars, " / " ) << "\".";
          continue;
        }
        if ( vars.size() == 3 && profile ) {
          const auto title = utils::format( "(%s / %s / %s) correlation;%s;%s;%s;d^{3}#sigma/d(%s)/d(%s)/d(%s) (pb/bin)",
            vars[0].c_str(), vars[1].c_str(), vars[2].c_str(),
            vars[0].c_str(), vars[1].c_str(), vars[2].c_str(),
            vars[0].c_str(), vars[1].c_str(), vars[2].c_str() );
          profiles2d_.emplace_back( std::make_pair( vars, YODA::Profile2D( nbins_x, min_x, max_x, nbins_y, min_y, max_y, key, title ) ) );
          CG_INFO( "YODAHistsHandler" )
            << "Booking a 2D profile"
            << " with " << utils::s( "bin", nbins_x+nbins_y )
            << " between (" << min_x << ", " << min_y << ")"
            << " and (" << max_x << ", " << max_y << ")"
            << " for \"" << utils::merge( vars, " / " ) << "\".";
          continue;
        }
      }
    }

    template<typename T>
    YODAHistsHandler<T>::~YODAHistsHandler()
    {
      std::vector<const YODA::AnalysisObject*> obj;
      //--- finalisation of the output file
      for ( const auto& hist : hists1d_ )
        obj.emplace_back( &hist.second );
      for ( const auto& hist : hists2d_ )
        obj.emplace_back( &hist.second );
      for ( const auto& hist : profiles1d_ )
        obj.emplace_back( &hist.second );
      for ( const auto& hist : profiles2d_ )
        obj.emplace_back( &hist.second );
      obj.emplace_back( &weight_cnt_ );
      T::write( file_, obj );
    }

    template<typename T> void
    YODAHistsHandler<T>::operator<<( const Event& ev )
    {
      //--- increment the corresponding histograms
      for ( auto& h_var : hists1d_ )
        h_var.second.fillBin( browser_.get( ev, h_var.first ), xsec_ );
      for ( auto& h_var : hists2d_ )
        h_var.second.fillBin(
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ), xsec_ );
      for ( auto& h_var : profiles1d_ )
        h_var.second.fill(
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ), xsec_ );
      for ( auto& h_var : profiles2d_ )
        h_var.second.fill(
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ),
          browser_.get( ev, h_var.first[2] ), xsec_ );
      weight_cnt_.fill( ev.weight );
    }
  }
}

typedef cepgen::io::YODAHistsHandler<YODA::WriterYODA> YodaOutputHandler;
typedef cepgen::io::YODAHistsHandler<YODA::WriterAIDA> YodaAidaOutputHandler;
typedef cepgen::io::YODAHistsHandler<YODA::WriterFLAT> YodaFlatOutputHandler;
REGISTER_IO_MODULE( "yoda", YodaOutputHandler )
REGISTER_IO_MODULE( "yoda_aida", YodaAidaOutputHandler )
REGISTER_IO_MODULE( "yoda_flat", YodaFlatOutputHandler )
