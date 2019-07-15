#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
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
        short extractVariableProperties( const std::string& );
        /// Retrieve a named variable from a particle
        double variable( const Particle&, const std::string& ) const;
        /// Retrieve a named variable from the whole event
        double variable( const Event&, const std::string& ) const;

        static const std::regex rgx_select_id_, rgx_select_role_;
        static constexpr double INVALID_OUTPUT = -999.;
        static constexpr size_t PLOT_WIDTH = 50;

        std::ofstream file_;
        const std::vector<std::string> variables_;
        const bool print_banner_, print_variables_;
        const ParametersList hist_variables_;
        const std::string separator_;

        //--- variables definition
        std::unordered_map<short,std::string> variables_name_;
        std::unordered_map<short,bool> variable_stored_;

        typedef std::pair<unsigned short,std::string> IndexedVariable;
        std::unordered_map<short,std::vector<IndexedVariable> > variables_per_id_;
        std::unordered_map<Particle::Role,std::vector<IndexedVariable> > variables_per_role_;
        std::vector<IndexedVariable> variables_for_event_;
        unsigned short num_vars_;
        std::ostringstream oss_vars_;

        double xsec_;

        //--- auxiliary helper maps
        const std::unordered_map<std::string,Particle::Role> role_str_ = {
          { "ib1", Particle::Role::IncomingBeam1 }, { "ib2", Particle::Role::IncomingBeam2 },
          { "ob1", Particle::Role::OutgoingBeam1 }, { "ob2", Particle::Role::OutgoingBeam2 },
          { "pa1", Particle::Role::Parton1 }, { "pa2", Particle::Role::Parton2 },
          { "cs",  Particle::Role::CentralSystem },
          { "int", Particle::Role::Intermediate }
        };
        typedef double( Particle::Momentum::*pMethod )(void) const;
        /// Mapping of string variables to momentum getter methods
        const std::unordered_map<std::string,pMethod> m_mom_str_ = {
          { "px",  &Particle::Momentum::px },
          { "py",  &Particle::Momentum::py },
          { "pz",  &Particle::Momentum::pz },
          { "pt",  &Particle::Momentum::pt },
          { "eta", &Particle::Momentum::eta },
          { "phi", &Particle::Momentum::phi },
          { "m",   &Particle::Momentum::mass },
          { "e",   &Particle::Momentum::energy },
          { "p",   &Particle::Momentum::p },
          { "pt2", &Particle::Momentum::pt2 },
          { "th",  &Particle::Momentum::theta },
          { "y",   &Particle::Momentum::rapidity }
        };

        //--- kinematic variables
        double sqrts_;
        unsigned long num_evts_;
        struct gsl_histogram_deleter
        {
          void operator()( gsl_histogram* h ) {
            gsl_histogram_free( h );
          }
        };
        std::unordered_map<short,std::unique_ptr<gsl_histogram,gsl_histogram_deleter> > hists_;
    };

    const std::regex TextHandler::rgx_select_id_( "(\\w+)\\((\\d+)\\)" );
    const std::regex TextHandler::rgx_select_role_( "(\\w+)\\(([a-z]+\\d?)\\)" );

    TextHandler::TextHandler( const ParametersList& params ) :
      GenericExportHandler( "text" ),
      file_           ( params.get<std::string>( "filename", "output.txt" ) ),
      variables_      ( params.get<std::vector<std::string> >( "variables" ) ),
      print_banner_   ( params.get<bool>( "saveBanner", true ) ),
      print_variables_( params.get<bool>( "saveVariables", true ) ),
      hist_variables_ ( params.get<ParametersList>( "histVariables" ) ),
      separator_      ( params.get<std::string>( "separator", "\t" ) ),
      num_vars_( 0 ), xsec_( 1. )
    {
      oss_vars_.clear();
      std::string sep;
      for ( const auto& var : variables_ ) {
        auto id = extractVariableProperties( var );
        if ( id >= 0 ) {
          oss_vars_ << sep << var, sep = separator_;
          variable_stored_[id] = true;
        }
      }
      for ( const auto& var : hist_variables_.keys() ) {
        auto id = extractVariableProperties( var );
        if ( id < 0 )
          continue;
        const auto& hvar = hist_variables_.get<ParametersList>( var );
        const int nbins = hvar.get<int>( "nbins", 10 );
        const double min = hvar.get<double>( "low", 0. ), max = hvar.get<double>( "high", 1. );
        hists_[id].reset( gsl_histogram_alloc( nbins ) );
        gsl_histogram_set_ranges_uniform( hists_[id].get(), min, max );
        CG_INFO( "TextHandler" )
          << "Booking a histogram with " << nbins << " bin" << utils::s( nbins )
          << " between " << min << " and " << max << " for \"" << var << "\".";
      }
    }

    TextHandler::~TextHandler()
    {
      for ( const auto& hs : hists_ ) {
        const auto& hist = hs.second.get();
        gsl_histogram_scale( hist, xsec_/( num_evts_+1 ) );
        const double max_bin = gsl_histogram_max_val( hist );
        const double inv_max_bin = max_bin > 0. ? 1./max_bin : 0.;
        CG_INFO( "TextHandler" )
          << "plot of \"" << variables_name_.at( hs.first ) << "\"\n"
          << std::string( 11, ' ' )
          << Form( "%-5.2f", gsl_histogram_min_val( hist ) )
          << std::string( PLOT_WIDTH-12, ' ' )
          << Form( "%5.2f", gsl_histogram_max_val( hist ) ) << " pb\n"
          << std::string( 11, ' ' )
          << std::string( PLOT_WIDTH+1, '.' );
        for ( size_t i = 0; i < gsl_histogram_bins( hist ); ++i ) {
          double min, max;
          gsl_histogram_get_range( hist, i, &min, &max );
          const int val = gsl_histogram_get( hist, i )*PLOT_WIDTH*inv_max_bin;
          CG_LOG( "TextHandler" ) << "["
            << std::setw( 4 ) << std::setprecision( 4 ) << min << ","
            << std::setw( 4 ) << std::setprecision( 4 ) << max << "):"
            << std::string( val, '*' );
        }
        //gsl_histogram_fprintf( stdout, hist, "%g", "%g" );
      }
      file_.close();
    }

    void
    TextHandler::initialise( const Parameters& params )
    {
      sqrts_ = params.kinematics.sqrtS();
      num_evts_ = 0ul;
      if ( print_banner_ )
        file_ << banner( params, "#" ) << "\n";
      if ( print_variables_ )
        file_ << "# " << oss_vars_.str() << "\n";
    }

    void
    TextHandler::operator<<( const Event& ev )
    {
      std::vector<double> vars( num_vars_ );
      //--- extract and order the variables to be retrieved
      //--- particle-level variables (indexed by integer id)
      for ( const auto& id_vars : variables_per_id_ ) {
        const auto& part = ev[id_vars.first];
        //--- loop over the list of variables for this particle
        for ( const auto& var : id_vars.second )
          vars[var.first] = variable( part, var.second );
      }
      //--- particle-level variables (indexed by role)
      for ( const auto& role_vars : variables_per_role_ ) {
        const auto& part = ev[role_vars.first][0];
        //--- loop over the list of variables for this particle
        for ( const auto& var : role_vars.second )
          vars[var.first] = variable( part, var.second );
      }
      //--- event-level variables
      for ( const auto& var : variables_for_event_ )
        vars[var.first] = variable( ev, var.second );
      //--- write down the variables list in the file
      std::string sep;
      unsigned short i = 0;
      for ( const auto& var : vars ) {
        if ( variable_stored_.count( i ) > 0 && variable_stored_.at( i ) )
          file_ << sep << var, sep = separator_;
        if ( hists_.count( i ) > 0 )
          gsl_histogram_increment( hists_.at( i ).get(), var );
        ++i;
      }
      file_ << "\n";
      ++num_evts_;
    }

    double
    TextHandler::variable( const Particle& part, const std::string& var ) const
    {
      if ( m_mom_str_.count( var ) ) {
        auto meth = m_mom_str_.at( var );
        return ( part.momentum().*meth )();
      }
      if ( var == "xi"  ) return 1.-part.energy()*2./sqrts_;
      if ( var == "pdg" ) return (double)part.integerPdgId();
      if ( var == "charge" ) return part.charge();
      if ( var == "status" ) return (double)part.status();
      CG_WARNING( "TextHandler" )
        << "Failed to retrieve variable \"" << var << "\".";
      return INVALID_OUTPUT;
    }

    double
    TextHandler::variable( const Event& ev, const std::string& var ) const
    {
      if ( var == "np" )
        return (double)ev.size();
      if ( var == "nev" )
        return (double)num_evts_+1;
      if ( var == "nob1" || var == "nob2" ) {
        unsigned short out = 0.;
        for ( const auto& part : ev[
          var == "nob1"
          ? Particle::Role::OutgoingBeam1
          : Particle::Role::OutgoingBeam2
        ] )
          if ( (int)part.status() > 0 )
            out++;
        return (double)out;
      }
      if ( var == "tgen" )
        return ev.time_generation;
      if ( var == "ttot" )
        return ev.time_total;
      CG_WARNING( "TextHandler" )
        << "Failed to retrieve the event-level variable \"" << var << "\".";
      return INVALID_OUTPUT;
    }

    short
    TextHandler::extractVariableProperties( const std::string& var )
    {
      const auto& vn = std::find_if( variables_name_.begin(), variables_name_.end(),
        [&var]( auto&& p ) { return p.second == var; } );
      if ( vn != variables_name_.end() )
        return vn->first;
      std::smatch sm;
      if ( std::regex_match( var, sm, rgx_select_id_ ) )
        variables_per_id_[std::stod( sm[2].str() )].emplace_back( std::make_pair( num_vars_, sm[1].str() ) );
      else if ( std::regex_match( var, sm, rgx_select_role_ ) ) {
        const auto& str_role = sm[2].str();
        if ( role_str_.count( str_role ) == 0 ) {
          CG_WARNING( "TextHandler" )
            << "Invalid particle role retrieved from configuration: \"" << str_role << "\".\n\t"
            << "Skipping the variable \"" << var << "\" in the output module.";
          return -1;
        }
        variables_per_role_[role_str_.at( str_role )].emplace_back( std::make_pair( num_vars_, sm[1].str() ) );
      }
      else // event-level variables
        variables_for_event_.emplace_back( std::make_pair( num_vars_, var ) );
      variables_name_[num_vars_] = var;
      return num_vars_++;
    }

  }
}

REGISTER_IO_MODULE( text, TextHandler )
