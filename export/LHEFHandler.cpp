#include "LHEFHandler.h"

#ifdef HEPMC_VERSION3

namespace OutputHandler
{
  LHEFHandler::LHEFHandler( const char* filename ) :
    HepMCHandler( filename, ExportHandler::LHE ),
    lhe_output_( std::make_unique<LHEF::Writer>( filename ) )
  {
    lhe_output_->init();
    HepMC::HEPEVT_Wrapper::set_hepevt_address( (char*)&hepevt_buf_ );
    HepMC::HEPEVT_Wrapper::zero_everything();
    lhe_output_->initComments() << "Sample created using CepGen.";
    lhe_output_->heprup.NPRUP = 1;
    lhe_output_->heprup.resize();
    lhe_output_->heprup.LPRUP[0] = -1;
    lhe_output_->init();
  }

  void
  LHEFHandler::operator<<( const Event* ev )
  {
    fillEvent( ev );
    if ( !event.get() or !HepMC::HEPEVT_Wrapper::GenEvent_to_HEPEVT( event.get() ) ) {
      throw Exception( __PRETTY_FUNCTION__, "Failed to retrieve the HepMC event to be stored!", FatalError );
    }
    LHEF::HEPEUP out;
    out.XWGTUP = 1.;
    out.XPDWUP = std::pair<double,double>( 0., 0. );
    out.SCALUP = 0.;
    out.AQEDUP = 0.;
    out.AQCDUP = 0.;
    out.NUP = hepevt_buf_.nhep;
    out.resize();
    for ( unsigned short ip=0; ip<hepevt_buf_.nhep; ip++ ) {
      out.IDUP[ip] = hepevt_buf_.idhep[ip]; // PDG id
      out.ISTUP[ip] = hepevt_buf_.isthep[ip]; // status code
      out.MOTHUP[ip] = std::pair<int,int>( hepevt_buf_.jmohep[ip][0], hepevt_buf_.jmohep[ip][1] ); // mothers
      out.ICOLUP[ip] = std::pair<int,int>( 0, 0 );
      out.PUP[ip] = std::vector<double>( hepevt_buf_.phep[ip][0], hepevt_buf_.phep[ip][4] ); // momentum
      std::cout << "before dumping!" << std::endl;
      for ( unsigned int i=0; i<5; i++ ) { std::cout << "pup[" << i << "]=" << out.PUP[ip][i] << std::endl;}
      std::cout << "finished dumping!" << std::endl;
      //out.VTIMUP.emplace_back( 5 ); // invariant lifetime
      //out.SPINUP.push_back( 0. );
    }
    out.print( std::cout );
    lhe_output_->hepeup = out;
    lhe_output_->writeEvent();
    event->clear();
  }
}

#endif
