#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Event/EventBrowser.h"
#include "CepGen/Parameters.h"

#include "CepGen/Version.h"

#include <gsl/gsl_histogram.h>

#include <iomanip>
#include <fstream>
#include <regex>

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the generic text file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class TextHandler : public GenericExportHandler
    {
      public:
        explicit TextHandler( const ParametersList& );
        ~TextHandler();

        void initialise( const Parameters& ) override;
        void setCrossSection( double xsec, double ) override { xsec_ = xsec; }
        void operator<<( const Event& ) override;

      private:
        std::string textHistogram( const std::string&, const gsl_histogram* ) const;

        static constexpr size_t PLOT_WIDTH = 50;
        static constexpr char PLOT_CHAR = '#';

        std::ofstream file_, hist_file_;
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
        struct gsl_histogram_deleter
        {
          void operator()( gsl_histogram* h ) {
            gsl_histogram_free( h );
          }
        };
        std::unordered_map<std::string,std::unique_ptr<gsl_histogram,gsl_histogram_deleter> > hists_;
    };

    TextHandler::TextHandler( const ParametersList& params ) :
      GenericExportHandler( "text" ),
      file_          ( params.get<std::string>( "filename", "output.txt" ) ),
      variables_     ( params.get<std::vector<std::string> >( "variables" ) ),
      save_banner_   ( params.get<bool>( "saveBanner", true ) ),
      save_variables_( params.get<bool>( "saveVariables", true ) ),
      show_hists_    ( params.get<bool>( "showHistograms", true ) ),
      save_hists_    ( params.get<bool>( "saveHistograms", false ) ),
      separator_     ( params.get<std::string>( "separator", "\t" ) ),
      xsec_( 1. )
    {
      //--- first extract list of variables to store in output file
      oss_vars_.clear();
      std::string sep;
      for ( const auto& var : variables_ )
        oss_vars_ << sep << var, sep = separator_;
      //--- then extract list of variables to be plotted in histogram
      const auto& hist_vars = params.get<ParametersList>( "histVariables" );
      for ( const auto& var : hist_vars.keys() ) {
        const auto& hvar = hist_vars.get<ParametersList>( var );
        const int nbins = hvar.get<int>( "nbins", 10 );
        const double min = hvar.get<double>( "low", 0. ), max = hvar.get<double>( "high", 1. );
        hists_[var].reset( gsl_histogram_alloc( nbins ) );
        gsl_histogram_set_ranges_uniform( hists_[var].get(), min, max );
        CG_INFO( "TextHandler" )
          << "Booking a histogram with " << nbins << " bin" << utils::s( nbins )
          << " between " << min << " and " << max << " for \"" << var << "\".";
      }
      if ( save_hists_ && !hists_.empty() )
        hist_file_.open( "lastrun.hists.txt" );
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
      ++num_evts_;
    }

    std::string
    TextHandler::textHistogram( const std::string& var, const gsl_histogram* hist ) const
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
        << sep << Form( "%-5.2f", gsl_histogram_min_val( hist ) )
        << std::string( PLOT_WIDTH-8, ' ' )
        << Form( "%5.2f", gsl_histogram_max_val( hist ) ) << "\n"
        << sep << std::string( PLOT_WIDTH+2, '.' ); // abscissa axis
      for ( size_t i = 0; i < nbins; ++i ) {
        double min, max;
        gsl_histogram_get_range( hist, i, &min, &max );
        const double value = gsl_histogram_get( hist, i );
        const int val = value*PLOT_WIDTH*inv_max_bin;
        os
          << "\n" << Form( "[%7.2f,%7.2f):", min, max )
          << std::string( val, PLOT_CHAR ) << std::string( PLOT_WIDTH-val, ' ' )
          << ": " << Form( "%6.2f", value );
      }
      const double bin_width = ( gsl_histogram_max( hist )-gsl_histogram_min( hist ) )/nbins;
      os
        << "\n"
        << Form( "%15s", var.c_str() ) << ":" << std::string( PLOT_WIDTH, '.' ) << ":\n" // 2nd abscissa axis
        << "\t("
        << "bin width=" << bin_width << " unit" << utils::s( (int)bin_width ) << ", "
        << "mean=" << gsl_histogram_mean( hist ) << ", "
        << "st.dev.=" << gsl_histogram_sigma( hist )
        << ")";
      return os.str();
    }
  }
}

REGISTER_IO_MODULE( text, TextHandler )
