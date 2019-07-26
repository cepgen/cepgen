#include "CepGen/IO/ExportHandler.h"

#include "CepGen/Core/Exception.h"
#include "CepGen/Core/ParametersList.h"
#include "CepGen/Core/utils.h"

#include "CepGen/Event/Event.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Parameters.h"

#include "promc/ProMCBook.h"

namespace cepgen
{
  namespace io
  {
    /**
     * \brief Handler for the generic ROOT file output
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
        std::unique_ptr<ProMCBook> file_;
        double xsec_, xsec_err_;
    };

    ProMCHandler::ProMCHandler( const ParametersList& params ) :
      GenericExportHandler( "promc" ),
      file_( new ProMCBook( params.get<std::string>( "filename", "output.promc" ).c_str(), "w" ) ),
      xsec_( -1. ), xsec_err_( -1. )
    {}

    ProMCHandler::~ProMCHandler()
    {
      ProMCStat stat;
      stat.set_cross_section_accumulated( xsec_ );
      stat.set_cross_section_error_accumulated( xsec_err_ );
      file_->setStatistics( stat );
      file_->close();
    }

    void
    ProMCHandler::initialise( const Parameters& params )
    {
      file_->setDescription( params.generation().maxgen, banner( params ) );
      ProMCHeader hdr;
      hdr.set_momentumunit( 1.e6 ); // in units of keV -> GeV
      hdr.set_lengthunit( 1.e3 ); //FIXME
      for ( const auto& pdg : PDG::get().particles() ) {
        auto data = hdr.add_particledata();
        data->set_id( (int)pdg );
        data->set_mass( PDG::get().mass( pdg ) );
        data->set_name( PDG::get().name( pdg ) );
        data->set_width( PDG::get().width( pdg ) );
        data->set_charge( PDG::get().charge( pdg ) );
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

      for ( const auto& par : ev.particles() ) {
        auto part = event.mutable_particles();
        part->add_id( par.integerPdgId() );
        part->add_status( (unsigned int)par.status() );
        //--- kinematics
        part->add_px( int( par.momentum().px()*1e6 ) );
        part->add_py( int( par.momentum().py()*1e6 ) );
        part->add_pz( int( par.momentum().pz()*1e6 ) );
        part->add_energy( int( par.energy()*1e6 ) );
        part->add_mass( int( par.mass()*1e6 ) );
        part->add_barcode( 0 );
        //--- parentage
        const auto& daugh = par.daughters(), &moth = par.mothers();
        part->add_daughter1( daugh.empty() ? 0 : *daugh.begin() );
        part->add_daughter2( daugh.size() > 1 ? *daugh.rbegin() : 0 );
        part->add_mother1( moth.empty() ? 0 : *moth.begin() );
        part->add_mother2( moth.size() > 1 ? *moth.rbegin() : 0 );
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

REGISTER_IO_MODULE( promc, ProMCHandler )
