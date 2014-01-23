#include "eventslist.h"

EventsList::EventsList(std::ofstream *of_, /*const HEPRUP *heprup_,*/ const int dumpevery_)
{
  std::cout << "[EventsList::EventsList] [DEBUG] Events list constructed" << std::endl;
  //_ev = new std::vector<Event>;
  this->_dump_every = dumpevery_;
  this->_lheof = of_;

  this->_lhe.str("");
  this->_lhe << "<LesHouchesEvents version=\"1.0\">" << std::endl;
  this->_lhe << "<header>This file was created from the output of the LPAIR++ generator</header>" << std::endl;
  /*this->_lhe << "<init>" << std::endl
	     << heprup_->idbmup[0] << " "
	     << heprup_->idbmup[1] << " "
	     << heprup_->ebmup[0] << " "
	     << heprup_->ebmup[1] << " "
	     << heprup_->pdfgup[0] << " "
	     << heprup_->pdfgup[1] << " "
	     << heprup_->pdfsup[0] << " "
	     << heprup_->pdfsup[1] << " "
	     << heprup_->idwtup << " "
	     << heprup_->nprup << std::endl;
  for (int i=0; i<heprup_->nprup; i++) {
    this->_lhe << heprup_->xsecup[i] << " "
	       << heprup_->xerrup[i] << " "
	       << heprup_->xmaxup[i] << " "
	       << heprup_->lprup[i] << std::endl;
  }
  //<< std::setprecision(2) << _ip.in1p << " "
  //<< std::setprecision(2) << _ip.in2p << " "
  //      << "0 0 10042 10042 2 1" << std::endl
  //	      << "0.10508723460E+01 0.96530000000E-02 0.26731120000E-03 0" << std::endl
  this->_lhe << "</init>" << std::endl;*/

  std::cout << this->_lhe.str() << std::endl;

}

EventsList::~EventsList()
{
  std::cout << "[EventsList::~EventsList] [DEBUG] Events list destroyed" << std::endl;

  this->_lhe << "</LesHouchesEvents>" << std::endl;
  this->DumpEvents();
  //delete _ev;
}

/*EventsList&
EventsList::operator=(const EventsList &ev_)
{
  this->_lhe = ev_._lhe;
  this->_ev = ev_._ev;
  return *this;
  }*/

void
EventsList::AddEvent(Event *ev_)
{
  //#ifdef DEBUG
  std::cout << "[EventsList::AddEvent] [DEBUG] New event added to the list (" << this->_ev.size()+1 << " elements)" << std::endl;
  //#endif
  this->_ev.push_back(*ev_);
  this->_lhe << ev_->GetLHERecord();
  if (this->NumEvents()!=0 and this->NumEvents()%this->_dump_every==0) {
    // dump into the output file
    this->DumpEvents();
    this->_lhe.str("");
  }
}

void
EventsList::DumpEvents()
{
  *this->_lheof << this->_lhe.str();
}

void
EventsList::Info()
{
  std::cout << "[EventsList::Info]" << std::endl
	    << this->NumEvents() << " events stored in the list"
	    << std::endl;
  std::vector<Event>::iterator it;
  for (it=this->_ev.begin(); it!=this->_ev.end(); it++) {
    (*it).Dump();
  }
}
