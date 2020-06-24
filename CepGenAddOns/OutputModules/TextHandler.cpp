#include "CepGen/Core/ExportModule.h"
#include "CepGen/Modules/ExportModuleFactory.h"

#include "CepGen/Core/Exception.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Parameters.h"

#include "CepGen/Utils/String.h"
#include "CepGen/Version.h"

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#include <iomanip>
#include <fstream>

namespace
{
  //--- 1D histograms
  struct gsl_histogram_deleter
  {
    void operator()( gsl_histogram* h ) {
      gsl_histogram_free( h );
    }
  };
  typedef std::unique_ptr<gsl_histogram,gsl_histogram_deleter> gsl_histogram_ptr;

  //--- 2D histograms
  struct gsl_histogram2d_deleter
  {
    void operator()( gsl_histogram2d* h ) {
      gsl_histogram2d_free( h );
    }
  };
  typedef std::unique_ptr<gsl_histogram2d,gsl_histogram2d_deleter> gsl_histogram2d_ptr;
}

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class TextHandler : public ExportModule
    {
      public:
        explicit TextHandler( const ParametersList& );
        ~TextHandler();

        void initialise( const Parameters& ) override;
        void setCrossSection( double xsec, double ) override { xsec_ = xsec; }
        void operator<<( const Event& ) override;

      private:
        static std::string textHistogram( const std::string&, const gsl_histogram* );
        static std::string textHistogram( const std::vector<std::string>&, const gsl_histogram2d* );

        static constexpr size_t PLOT_WIDTH = 50;
        static constexpr char PLOT_CHAR = '#';
        static const std::vector<char> PLOT_2D_CHARS;

        std::ofstream file_, hist_file_;
        std::string hist_filename_;
        //--- variables definition
        const std::vector<std::string> variables_;
        const bool save_banner_, save_variables_;
        const bool show_hists_, save_hists_;
        const std::string separator_;

        const utils::EventBrowser browser_;

        std::ostringstream oss_vars_;

        double xsec_;

        //--- kinematic variables
        double sqrts_;
        unsigned long num_evts_;
        std::vector<std::pair<std::string,gsl_histogram_ptr> > hists_;
        std::vector<std::pair<std::vector<std::string>,gsl_histogram2d_ptr> > hists2d_;
    };

    const std::vector<char>
    TextHandler::PLOT_2D_CHARS = {' ', '.', 'o', 'O', '0', '#'};

    TextHandler::TextHandler( const ParametersList& params ) :
      ExportModule( params ),
      file_          ( params.get<std::string>( "filename", "output.txt" ) ),
      hist_filename_ ( params.get<std::string>( "histFilename", "output.hists.txt" ) ),
      variables_     ( params.get<std::vector<std::string> >( "variables" ) ),
      save_banner_   ( params.get<bool>( "saveBanner", true ) ),
      save_variables_( params.get<bool>( "saveVariables", true ) ),
      show_hists_    ( params.get<bool>( "showHistograms", true ) ),
      save_hists_    ( params.get<bool>( "saveHistograms", false ) ),
      separator_     ( params.get<std::string>( "separator", "\t" ) ),
      xsec_( 1. ), sqrts_( 0. ), num_evts_( 0ul )
    {
      //--- first extract list of variables to store in output file
      oss_vars_.clear();
      std::string sep;
      for ( const auto& var : variables_ )
        oss_vars_ << sep << var, sep = separator_;
      //--- then extract list of variables to be plotted in histogram
      const auto& hist_vars = params.get<ParametersList>( "histVariables" );
      for ( const auto& key : hist_vars.keys() ) {
        const auto& vars = utils::split( key, ':' );
        if ( vars.size() < 1 || vars.size() > 2 )
          throw CG_FATAL( "ROOTHistsHandler" )
            << "Invalid number of variables to correlate for '" << key << "'!";

        const auto& hvar = hist_vars.get<ParametersList>( key );
        const int nbins_x = hvar.get<int>( "nbinsX", hvar.get<int>( "nbins", PLOT_WIDTH/2 ) );
        const double min_x = hvar.get<double>( "lowX", hvar.get<double>( "low", 0. ) );
        const double max_x = hvar.get<double>( "highX", hvar.get<double>( "high", 1. ) );
        if ( vars.size() == 1 ) { // 1D histogram
          auto hist = gsl_histogram_alloc( nbins_x );
          gsl_histogram_set_ranges_uniform( hist, min_x, max_x );
          hists_.emplace_back( std::make_pair( key, gsl_histogram_ptr( hist ) ) );
          CG_INFO( "TextHandler" )
            << "Booking a 1D histogram with " << utils::s( "bin", nbins_x )
            << " between " << min_x << " and " << max_x << " for \"" << key << "\".";
        }
        else if ( vars.size() == 2 ) { // 2D histogram
          const int nbins_y = hvar.get<int>( "nbinsY", PLOT_WIDTH );
          const double min_y = hvar.get<double>( "lowY", 0. );
          const double max_y = hvar.get<double>( "highY", 1. );
          auto hist = gsl_histogram2d_alloc( nbins_x, nbins_y );
          gsl_histogram2d_set_ranges_uniform( hist, min_x, max_x, min_y, max_y );
          hists2d_.emplace_back( std::make_pair( vars, gsl_histogram2d_ptr( hist ) ) );
          CG_INFO( "TextHandler" )
            << "Booking a 2D correlation plot with " << utils::s( "bin", nbins_x+nbins_y )
            << " between (" << min_x << ", " << max_x << ") and (" << min_y << ", " << max_y << ")"
            << " for \"" << key << "\".";
        }
      }
      if ( save_hists_ && !hists_.empty() )
        hist_file_.open( hist_filename_ );
    }

    TextHandler::~TextHandler()
    {
      //--- finalisation of the output file
      file_.close();
      //--- histograms printout
      if ( !show_hists_ && !save_hists_ )
        return;
      for ( const auto& h_var : hists_ ) {
        const auto& hist = h_var.second.get();
        gsl_histogram_scale( hist, xsec_/( num_evts_+1 ) );
        const auto& h_out = textHistogram( h_var.first, hist );
        if ( show_hists_ )
          CG_INFO( "TextHandler" ) << h_out;
        if ( save_hists_ )
          hist_file_ << "\n" << h_out << "\n";
      }
      for ( const auto& h_var : hists2d_ ) {
        const auto& hist = h_var.second.get();
        //gsl_histogram_scale( hist, xsec_/( num_evts_+1 ) );
        const auto& h_out = textHistogram( h_var.first, hist );
        if ( show_hists_ )
          CG_INFO( "TextHandler" ) << h_out;
        if ( save_hists_ )
          hist_file_ << "\n" << h_out << "\n";
      }
      if ( save_hists_ )
        CG_INFO( "TextHandler" )
          << "Saved " << utils::s( "histogram", hists_.size() )
          << " into \"" << hist_filename_ << "\".";
    }

    void
    TextHandler::initialise( const Parameters& params )
    {
      sqrts_ = params.kinematics.sqrtS();
      num_evts_ = 0ul;
      if ( save_banner_ )
        file_ << banner( params, "#" ) << "\n";
      if ( save_variables_ )
        file_ << "# " << oss_vars_.str() << "\n";
      if ( save_hists_ && !hists_.empty() )
        hist_file_ << banner( params, "#" ) << "\n";
    }

    void
    TextHandler::operator<<( const Event& ev )
    {
      //--- write down the variables list in the file
      std::string sep;
      for ( const auto& var : variables_ )
        file_ << sep << browser_.get( ev, var ), sep = separator_;
      file_ << "\n";
      //--- increment the corresponding histograms
      for ( const auto& h_var : hists_ )
        gsl_histogram_increment( h_var.second.get(), browser_.get( ev, h_var.first ) );
      for ( const auto& h_var : hists2d_ )
        gsl_histogram2d_increment( h_var.second.get(),
          browser_.get( ev, h_var.first[0] ),
          browser_.get( ev, h_var.first[1] ) );
      ++num_evts_;
    }

    std::string
    TextHandler::textHistogram( const std::string& var, const gsl_histogram* hist )
    {
      std::ostringstream os;
      const size_t nbins = gsl_histogram_bins( hist );
      const double max_bin = gsl_histogram_max_val( hist );
      const double inv_max_bin = max_bin > 0. ? 1./max_bin : 0.;
      const std::string sep( 17, ' ' );
      os
        << "plot of \"" << var << "\"\n"
        << sep << std::string( PLOT_WIDTH-15-var.size(), ' ' )
        << "d(sig)/d" << var << " (pb/bin)\n"
        << sep << utils::format( "%-5.2f", gsl_histogram_min_val( hist ) )
        << std::string( PLOT_WIDTH-11, ' ' )
        << utils::format( "%5.2e", gsl_histogram_max_val( hist ) ) << "\n"
        << sep << std::string( PLOT_WIDTH+2, '.' ); // abscissa axis
      for ( size_t i = 0; i < nbins; ++i ) {
        double min, max;
        gsl_histogram_get_range( hist, i, &min, &max );
        const double value = gsl_histogram_get( hist, i );
        const int val = value*PLOT_WIDTH*inv_max_bin;
        os
          << "\n" << utils::format( "[%7.2f,%7.2f):", min, max )
          << std::string( val, PLOT_CHAR ) << std::string( PLOT_WIDTH-val, ' ' )
          << ": " << utils::format( "%6.2e", value );
      }
      const double bin_width = ( gsl_histogram_max( hist )-gsl_histogram_min( hist ) )/nbins;
      os
        << "\n"
        << utils::format( "%17s", var.c_str() ) << ":" << std::string( PLOT_WIDTH, '.' ) << ":\n" // 2nd abscissa axis
        << "\t("
        << "bin width=" << bin_width << utils::s( " unit", (int)bin_width, false ) << ", "
        << "mean=" << gsl_histogram_mean( hist ) << ", "
        << "st.dev.=" << gsl_histogram_sigma( hist )
        << ")";
      return os.str();
    }

    std::string
    TextHandler::textHistogram( const std::vector<std::string>& vars, const gsl_histogram2d* hist )
    {
      std::ostringstream os;
      const size_t nbins_x = gsl_histogram2d_nx( hist );
      const size_t nbins_y = gsl_histogram2d_ny( hist );
      const double max_bin = gsl_histogram2d_max_val( hist );
      const double inv_max_bin = max_bin > 0. ? 1./max_bin : 0.;
      const std::string sep( 17, ' ' );
      const auto var = utils::merge( vars, "/" );
      os
        << "plot of \"" << var << "\"\n"
        << sep << std::string( (size_t)std::max( 0., nbins_y-15.-var.size() ), ' ' );
      os
        << "d^2(sig)/d" << vars.at( 0 ) << "/d" << vars.at( 1 ) << " (pb/bin)\n"
        << sep << utils::format( "%-5.2f", gsl_histogram2d_ymin( hist ) )
        << std::string( nbins_y-11, ' ' )
        << utils::format( "%5.2e", gsl_histogram2d_ymax( hist ) ) << "\n"
        << utils::format( "%17s", vars.at( 0 ).c_str() )
        << std::string( nbins_y+2, '.' ); // abscissa axis
      for ( size_t i = 0; i < nbins_x; ++i ) {
        double min_x, max_x;
        gsl_histogram2d_get_xrange( hist, i, &min_x, &max_x );
        os
          << "\n" << utils::format( "[%7.2f,%7.2f):", min_x, max_x );
        for ( size_t j = 0; j < nbins_y; ++j ) {
          const double value_norm = gsl_histogram2d_get( hist, i, j )*inv_max_bin;
          os << PLOT_2D_CHARS[ceil(value_norm*(PLOT_2D_CHARS.size()-1))];
        }
        os << ":";
      }
      std::vector<std::string> ylabels;
      for ( size_t j = 0; j < nbins_y; ++j ) {
        double min_y, max_y;
        gsl_histogram2d_get_yrange( hist, j, &min_y, &max_y );
        ylabels.emplace_back( utils::format( "%g", min_y ) );
      }
      struct stringlen {
        bool operator()( const std::string& a, const std::string& b ) {
          return a.size() < b.size();
        }
      };
      for ( size_t i = 0; i < std::max_element( ylabels.begin(), ylabels.end(), stringlen() )->size(); ++i ) {
        os << "\n" << sep << ":";
        for ( const auto& lab : ylabels )
          os << ( lab.size() > i ? lab.at( i ) : ' ' );
        os << ":";
      }
      const double bin_width_x = ( gsl_histogram2d_xmax( hist )-gsl_histogram2d_xmin( hist ) )/nbins_x;
      const double bin_width_y = ( gsl_histogram2d_ymax( hist )-gsl_histogram2d_ymin( hist ) )/nbins_y;
      os
        << "\n" << sep
        << ":" << std::string( nbins_y, '.' ) << ": " // 2nd abscissa axis
        << vars.at( 1 ) << "\n\t("
        << "x-axis: "
        << "bin width=" << bin_width_x << utils::s( " unit", (int)bin_width_x, false ) << ", "
        << "mean=" << gsl_histogram2d_xmean( hist ) << ","
        << "st.dev.=" << gsl_histogram2d_xsigma( hist ) << "\n\t"
        << " y-axis: "
        << "bin width=" << bin_width_y << utils::s( " unit", (int)bin_width_y, false ) << ", "
        << "mean=" << gsl_histogram2d_ymean( hist ) << ","
        << "st.dev.=" << gsl_histogram2d_ysigma( hist )
        << ")";
      return os.str();
    }
  }
}

REGISTER_IO_MODULE( "text", TextHandler )
