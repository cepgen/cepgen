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
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the generic ROOT file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class ROOTHistsHandler : public GenericExportHandler
    {
      public:
        explicit ROOTHistsHandler( const ParametersList& );
        ~ROOTHistsHandler();

        void initialise( const Parameters& ) override {}
        void setCrossSection( double xsec, double ) override { xsec_ = xsec; }
        void operator<<( const Event& ) override;

      private:
        std::unique_ptr<TFile> file_;
        std::vector<std::pair<std::string,TH1*> > hists1d_;
        std::vector<std::pair<std::vector<std::string>,TH2*> > hists2d_;
        std::vector<std::pair<std::vector<std::string>,TH3*> > hists3d_;
        std::vector<std::pair<std::vector<std::string>,TProfile*> > profiles1d_;
        std::vector<std::pair<std::vector<std::string>,TProfile2D*> > profiles2d_;

        const ParametersList variables_;

        double xsec_;
        const utils::EventBrowser browser_;
    };

    ROOTHistsHandler::ROOTHistsHandler( const ParametersList& params ) :
      GenericExportHandler( "root-hist" ),
      file_( TFile::Open( params.get<std::string>( "filename", "output.root" ).c_str(), "recreate" ) ),
      variables_( params.get<ParametersList>( "variables" ) ),
      xsec_( 1. )
    {
      //--- extract list of variables/correlations to be plotted in histograms
      for ( const auto& key : variables_.keys() ) {
        const auto& vars = split( key, ':' );
        if ( vars.size() < 1 || vars.size() > 3 )
          throw CG_FATAL( "ROOTHistsHandler" )
            << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvars = variables_.get<ParametersList>( key );
        int nbins_x = hvars.get<int>( "nbinsX", 10 );
        nbins_x = hvars.get<int>( "nbins", nbins_x );
        double min_x = hvars.get<double>( "lowX", 0. ), max_x = hvars.get<double>( "highX", 1. );
        min_x = hvars.get<double>( "low", min_x ), max_x = hvars.get<double>( "high", max_x );
        const bool profile = hvars.get<bool>( "profile", false );
        if ( vars.size() == 1 ) { // 1D histogram
          const auto title = Form( "%s;%s;d#sigma/d(%s) (pb/bin)", key.c_str(), key.c_str(), key.c_str() );
          hists1d_.emplace_back( std::make_pair( key,
            new TH1D( key.c_str(), title.c_str(), nbins_x, min_x, max_x ) ) );
          CG_INFO( "ROOTHistsHandler" )
            << "Booking a 1D histogram with " << utils::s( "bin", nbins_x )
            << " between " << min_x << " and " << max_x
            << " for \"" << key << "\".";
          continue;
        }
        const int nbins_y = hvars.get<int>( "nbinsY", 10 );
        const double min_y = hvars.get<double>( "lowY", 0. ), max_y = hvars.get<double>( "highY", 1. );
        if ( vars.size() == 2 ) { // 2D histogram / 1D profile
          const auto title = Form( "(%s / %s) correlation;%s;%s;d^{2}#sigma/d(%s)/d(%s) (pb/bin)",
            vars[0].c_str(), vars[1].c_str(),
            vars[0].c_str(), vars[1].c_str(),
            vars[0].c_str(), vars[1].c_str() );
          if ( profile )
            profiles1d_.emplace_back( std::make_pair( vars,
              new TProfile( key.c_str(), title.c_str(), nbins_x, min_x, max_x ) ) );
          else
            hists2d_.emplace_back( std::make_pair( vars,
              new TH2D( key.c_str(), title.c_str(), nbins_x, min_x, max_x, nbins_y, min_y, max_y ) ) );
          CG_INFO( "ROOTHistsHandler" )
            << "Booking a "
            << ( profile ? "1D profile" : "2D correlation plot" )
            << " with " << utils::s( "bin", nbins_x+nbins_y )
            << " between (" << min_x << ", " << min_y << ") and (" << max_x << ", " << max_y << ")"
            << " for \"" << merge( vars, " / " ) << "\".";
          continue;
        }
        const int nbins_z = hvars.get<int>( "nbinsZ", 10 );
        const double min_z = hvars.get<double>( "lowZ", 0. ), max_z = hvars.get<double>( "highZ", 1. );
        if ( vars.size() == 3 ) { // 3D histogram
          const auto title = Form( "(%s / %s / %s) correlation;%s;%s;%s;d^{3}#sigma/d(%s)/d(%s)/d(%s) (pb/bin)",
            vars[0].c_str(), vars[1].c_str(), vars[2].c_str(),
            vars[0].c_str(), vars[1].c_str(), vars[2].c_str(),
            vars[0].c_str(), vars[1].c_str(), vars[2].c_str() );
          if ( profile )
            profiles2d_.emplace_back( std::make_pair( vars,
              new TProfile2D( key.c_str(), title.c_str(), nbins_x, min_x, max_x, nbins_y, min_y, max_y ) ) );
          else
            hists3d_.emplace_back( std::make_pair( vars,
              new TH3D( key.c_str(), title.c_str(), nbins_x, min_x, max_x, nbins_y, min_y, max_y, nbins_z, min_z, max_z ) ) );
          CG_INFO( "ROOTHistsHandler" )
            << "Booking a "
            << ( profile ? "2D profile" : "3D correlation plot" )
            << " with " << utils::s( "bin", nbins_x+nbins_y+nbins_z )
            << " between (" << min_x << ", " << min_y << ", " << min_z << ")"
            << " and (" << max_x << ", " << max_y << ", " << max_z << ") "
            << "for \"" << merge( vars, " / " ) << "\".";
          continue;
        }
      }
    }

    ROOTHistsHandler::~ROOTHistsHandler()
    {
      //--- finalisation of the output file
      for ( const auto& hist : hists1d_ )
        hist.second->Write( hist.first.c_str() );
      for ( const auto& hist : hists2d_ )
        hist.second->Write( merge( hist.first, "_vs_" ).c_str() );
      for ( const auto& hist : hists3d_ )
        hist.second->Write( merge( hist.first, "_vs_" ).c_str() );
      for ( const auto& hist : profiles1d_ )
        hist.second->Write( merge( hist.first, "_vs_" ).c_str() );
      for ( const auto& hist : profiles2d_ )
        hist.second->Write( merge( hist.first, "_vs_" ).c_str() );
      // ROOT and its sumptuous memory management disallows the "delete" here
      file_->Close();
    }

    void
    ROOTHistsHandler::operator<<( const Event& ev )
    {
      //--- increment the corresponding histograms
      for ( const auto& h_var : hists1d_ )
        h_var.second->Fill(
          browser_.get( ev, h_var.first ), xsec_ );
      for ( const auto& h_var : hists2d_ )
        h_var.second->Fill(
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ), xsec_ );
      for ( const auto& h_var : hists3d_ )
        h_var.second->Fill(
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ),
          browser_.get( ev, h_var.first[2] ), xsec_ );
      for ( const auto& h_var : profiles1d_ )
        h_var.second->Fill(
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ), xsec_ );
      for ( const auto& h_var : profiles2d_ )
        h_var.second->Fill(
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ),
          browser_.get( ev, h_var.first[2] ), xsec_ );
    }
  }
}

REGISTER_IO_MODULE( "root_hist", ROOTHistsHandler )
