#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Parameters.h"
#include "CepGen/Version.h"

#include "ProMCBook.h"

#include <stdio.h>

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the ProMC file output
     * \author Laurent Forthomme <laurent.forthomme@cern.ch>
     * \date Jul 2019
     */
    class ProMCHandler : public GenericExportHandler
    {
      public:
        explicit ProMCHandler( const ParametersList& );
        ~ProMCHandler();

        void initialise( const Parameters& ) override;
        void setCrossSection( double xsec, double err ) override { xsec_ = xsec, xsec_err_ = err; }
        void operator<<( const Event& ) override;

      private:
        static constexpr double GEV_UNIT = 1.e6; // base unit in GEV_UNIT^-1 GeV = keV
        static constexpr double M_UNIT = 1.e3; // base unit in M^-1 m = mm
        static int inGeV( double val ) { return int( val*GEV_UNIT ); }

        std::unique_ptr<ProMCBook> file_;
        const bool compress_evt_;
        std::ofstream log_file_;
        double xsec_, xsec_err_;
    };

    ProMCHandler::ProMCHandler( const ParametersList& params ) :
      GenericExportHandler( "promc" ),
      file_( new ProMCBook( params.get<std::string>( "filename", "output.promc" ).c_str(), "w" ) ),
      compress_evt_( params.get<bool>( "compress", false ) ),
      log_file_( "logfile.txt" ),
      xsec_( -1. ), xsec_err_( -1. )
    {}

    ProMCHandler::~ProMCHandler()
    {
      ProMCStat stat;
      stat.set_cross_section_accumulated( xsec_ );
      stat.set_cross_section_error_accumulated( xsec_err_ );
      stat.set_luminosity_accumulated( event_num_/xsec_ );
      stat.set_ntried( event_num_ );
      stat.set_nselected( event_num_ );
      stat.set_naccepted( event_num_ );
      file_->setStatistics( stat );
      file_->close();
      //--- delete the log file once attached
      remove( "logfile.txt" );
    }

    void
    ProMCHandler::initialise( const Parameters& params )
    {
      file_->setDescription( params.generation().maxgen, "Sample generated using CepGen v"+version() );
      log_file_ << banner( params ) << "\n";
      ProMCHeader hdr;
      hdr.set_momentumunit( GEV_UNIT );
      hdr.set_lengthunit( M_UNIT ); // unused as for now
      for ( const auto& pdg : PDG::get().particles() ) {
        auto data = hdr.add_particledata();
        const auto& desc = PDG::get()( pdg );
        data->set_id( (int)pdg );
        data->set_mass( desc.mass );
        data->set_name( desc.name );
        data->set_width( desc.width );
        data->set_charge( desc.charge*1./3. );
      }
      hdr.set_id1( params.kinematics.incoming_beams.first.pdg );
      hdr.set_id2( params.kinematics.incoming_beams.second.pdg );
      hdr.set_pdf1( 0 );
      hdr.set_pdf2( 0 );
      hdr.set_x1( 0 );
      hdr.set_x2( 0 );
      hdr.set_ecm( params.kinematics.sqrtS() );
      file_->setHeader( hdr );
    }

    void
    ProMCHandler::operator<<( const Event& ev )
    {
      ProMCEvent event;
      auto evt = event.mutable_event();
      evt->set_number( event_num_++ );
      evt->set_process_id( 0 );
      evt->set_scale( ev[Particle::Role::Intermediate][0].mass() );
      evt->set_alpha_qed( constants::ALPHA_EM );
      evt->set_alpha_qcd( constants::ALPHA_QCD );
      evt->set_weight( 1. );

      unsigned short i = 0;
      const auto& parts = compress_evt_
        ? ev.compressed().particles()
        : ev.particles();
      for ( const auto& par : parts ) {
        auto part = event.mutable_particles();
        part->add_id( i++ );
        part->add_pdg_id( par.integerPdgId() );
        part->add_status( (unsigned int)par.status() );
        //--- kinematics
        part->add_px( inGeV( par.momentum().px() ) );
        part->add_py( inGeV( par.momentum().py() ) );
        part->add_pz( inGeV( par.momentum().pz() ) );
        part->add_energy( inGeV( par.energy() ) );
        part->add_mass( inGeV( par.mass() ) );
        part->add_barcode( 0 );
        //--- parentage
        const auto& daugh = par.daughters(), &moth = par.mothers();
        part->add_daughter1( daugh.empty() ? 0 : *daugh.begin()+1 );
        part->add_daughter2( daugh.size() > 1 ? *daugh.rbegin()+1 : 0 );
        part->add_mother1( moth.empty() ? 0 : *moth.begin()+1 );
        part->add_mother2( moth.size() > 1 ? *moth.rbegin()+1 : 0 );
        //--- vertex
        part->add_x( 0 );
        part->add_y( 0 );
        part->add_z( 0 );
        part->add_t( 0 );
      }
      file_->write( event );
    }
  }
}

REGISTER_IO_MODULE( "promc", ProMCHandler )
